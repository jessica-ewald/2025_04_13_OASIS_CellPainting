import polars as pl
from pathlib import Path

input_dir = Path("../01_snakemake/inputs")

dat = pl.read_parquet(input_dir / Path("profiles/cellprofiler/raw.parquet"))

dat = dat.filter(pl.col("Metadata_well_type") == "OASIS_cmpd_exposure")
cmpds = dat.select(["Metadata_Compound", "Metadata_Source", "Metadata_Plate", "Metadata_Well"])

print("Number of unique batches, plates, and well positions per compound:")
print(cmpds.group_by("Metadata_Compound").agg([pl.col("Metadata_Source").unique().len(), pl.col("Metadata_Plate").unique().len(), pl.col("Metadata_Well").unique().len()]))