import logging  # noqa: CPY001, D100

import polars as pl

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def call_hits(
    mtt_path: str,
    ldh_path: str,
    cc_path: str,
    hits_path: str,
) -> None:
    """Binarize assay curves."""
    # 1. Read in data
    mtt = pl.read_parquet(mtt_path)
    ldh = pl.read_parquet(ldh_path)
    cc = pl.read_parquet(cc_path)

    # 2. Get hits
    hits = mtt.select("Metadata_OASIS_ID")

    mtt_hits = mtt.filter(pl.col("all.pass") == True).select("Metadata_OASIS_ID").to_series().to_list()
    ldh_hits = ldh.filter(pl.col("all.pass") == True).select("Metadata_OASIS_ID").to_series().to_list()
    cc_hits = cc.filter(pl.col("all.pass") == True).select("Metadata_OASIS_ID").to_series().to_list()

    hits = hits.with_columns(
        pl.when(pl.col("Metadata_OASIS_ID").is_in(mtt_hits)).then(pl.lit(1)).otherwise(pl.lit(0)).alias("MTT"),
        pl.when(pl.col("Metadata_OASIS_ID").is_in(ldh_hits)).then(pl.lit(1)).otherwise(pl.lit(0)).alias("LDH"),
        pl.when(pl.col("Metadata_OASIS_ID").is_in(cc_hits)).then(pl.lit(1)).otherwise(pl.lit(0)).alias("cell_count"),
    )

    hits = hits.rename({"Metadata_OASIS_ID": "OASIS_ID"})

    # 3. Write out results
    hits.write_parquet(hits_path)
