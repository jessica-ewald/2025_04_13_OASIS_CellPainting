import os

import polars as pl
from tqdm import tqdm


def main() -> None:
    """Format CellProfiler profiles.

    Merge profiles from each plate into one file and add metadata.

    """
    input_profile_path = "../1_snakemake/inputs/profiles/cellprofiler/plates"
    meta_path = "../1_snakemake/inputs/metadata/metadata.parquet"
    output_profile_path = "../1_snakemake/inputs/profiles/cellprofiler/raw.parquet"

    meta = pl.read_parquet(meta_path)

    profiles = []
    plates = os.listdir(input_profile_path)
    plates = [i for i in plates if "plate_" in i]

    # Get column schema
    schema = pl.read_csv(f"{input_profile_path}/{plates[0]}", infer_schema_length=10000)
    meta_cols = [col for col in schema.columns if "Metadata" in col]
    schema = schema.with_columns([pl.col(col).cast(pl.Float64) for col in schema.columns if col not in meta_cols])
    schema = schema.schema

    # Read in data for each plate
    for plate in tqdm(plates):
        prof_path = f"{input_profile_path}/{plate}"
        profile = pl.read_csv(prof_path, schema=schema)
        profiles.append(profile)

    # Concat together
    data = pl.concat(profiles, how="vertical_relaxed")
    data = meta.join(data, on=["Metadata_Plate", "Metadata_Well"])
    data.write_parquet(output_profile_path)


if __name__ == "__main__":
    main()
