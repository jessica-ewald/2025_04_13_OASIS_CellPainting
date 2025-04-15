import os

from pathlib import Path
import numpy as np
import polars as pl

meta_keep = [
    "plate_map_name",
    "well_position",
    "BatchId",
    "Concentration(uM)"
]


def process_meta(input_meta_path: str, meta_nms: list) -> pl.DataFrame:
    """Process metadata.

    Select columns of interest and rename.



    """
    meta = pl.read_csv(input_meta_path)
    missing_cols = [i for i in meta_keep if i not in meta.columns]
    for mc in missing_cols:
        meta = meta.with_columns(
            pl.lit(None).alias(mc),
        )
    meta = meta.select(meta_keep)
    meta.columns = meta_nms

    return meta


def main() -> None:
    """Format metadata.

    Merge metadata from each plate into one file.

    """
    # Process metadata
    meta_nms = [f"Metadata_{i}" for i in meta_keep]
    meta_path = Path("../01_snakemake/inputs/metadata/plates/")
    meta = []

    plates = os.listdir(meta_path)
    for plate in plates:
        plate_path = meta_path / plate
        meta.append(process_meta(plate_path, meta_nms))

    meta = pl.concat(meta, how="vertical_relaxed")
    meta = meta.rename({
        "Metadata_plate_map_name": "Metadata_Plate",
        "Metadata_well_position": "Metadata_Well",
        "Metadata_BatchId": "Metadata_Broad_ID",
        "Metadata_Concentration(uM)": "Metadata_Concentration",
    }).with_columns(
        pl.when(pl.col("Metadata_Concentration").is_null()).then(pl.lit(0)).otherwise(pl.col("Metadata_Concentration")).alias("Metadata_Concentration")
    ).with_columns(
        pl.concat_str(["Metadata_Broad_ID", "Metadata_Concentration"], separator="_").alias("Metadata_Perturbation"),
    )

    compounds = meta.filter(pl.col("Metadata_Concentration") != 0).select("Metadata_Broad_ID").to_series().unique().to_list()

    meta_log10 = meta.filter(pl.col("Metadata_Concentration") == 0).with_columns(
        pl.lit(0).cast(pl.Float64).alias("Metadata_Log10Conc"),
    )

    # Convert non-zero concentrations to log-scale and shift to above zero
    for compound in compounds:
        temp = meta.filter(pl.col("Metadata_Broad_ID") == compound)
        concs = temp.select(pl.col("Metadata_Concentration")).to_series().sort().unique().to_list()

        if len(concs) > 1:
            shift_val = np.abs(np.log10(concs[1] / concs[0]))

            temp = temp.with_columns(
                (pl.col("Metadata_Concentration").log10() + shift_val).alias("Metadata_Log10Conc"),
            )
        else:
            temp = temp.with_columns(
                (pl.lit(None)).alias("Metadata_Log10Conc")
            )

        meta_log10 = pl.concat([meta_log10, temp], how="vertical")

    # need to merge with OASIS metadata
    id_map = pl.read_csv("../01_snakemake/inputs/annotations/OASIS_BRDID_map.csv").select([
        "OASIS_ID", "BROAD_ID"
    ]).rename({"OASIS_ID": "Metadata_OASIS_ID", "BROAD_ID": "Metadata_Broad_ID"})

    oasis = pl.read_csv("../01_snakemake/inputs/annotations/v5_oasis_03Sept2024_simple.csv").select([
        "OASIS_ID", "DTXSID", "PREFERRED_NAME", "INCHIKEY", "Purchased_Axiom_Medchemxpress"
    ]).rename({
        "OASIS_ID": "Metadata_OASIS_ID",
        "DTXSID": "Metadata_DTXSID",
        "PREFERRED_NAME": "Metadata_Compound",
        "INCHIKEY": "Metadata_Inchikey",
        "Purchased_Axiom_Medchemxpress": "Metadata_Axiom"
    }).filter(~pl.col("Metadata_OASIS_ID").is_null()).unique()

    meta_log10 = meta_log10.join(id_map, on="Metadata_Broad_ID", how="left")
    meta_log10 = meta_log10.join(oasis, on="Metadata_OASIS_ID", how="left")

    meta_log10 = meta_log10.with_columns(
        pl.when(pl.col("Metadata_Concentration") == 0).then(pl.lit("DMSO")).when((pl.col("Metadata_Concentration") > 0) & (pl.col("Metadata_Compound").is_null())).then(pl.lit("JUMP_control")).otherwise(pl.lit("OASIS_cmpd_exposure")).alias("Metadata_well_type")
    )

    # rename blank compound as DMSO or BRD ID
    meta_log10 = meta_log10.with_columns(
        pl.when(pl.col("Metadata_Concentration") == 0).then(pl.lit("DMSO")).otherwise(pl.col("Metadata_Compound")).alias("Metadata_Compound")
    )

    meta_log10 = meta_log10.with_columns(
        pl.when(pl.col("Metadata_Compound").is_null()).then(pl.col("Metadata_Broad_ID")).otherwise(pl.col("Metadata_Compound")).alias("Metadata_Compound")
    )

    meta_log10.write_parquet("../01_snakemake/inputs/metadata/metadata.parquet")

if __name__ == "__main__":
    main()
