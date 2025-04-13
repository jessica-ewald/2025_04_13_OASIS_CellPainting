import logging  # noqa: CPY001, D100

import polars as pl
import pycytominer

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def aggregate_compound(method: str, dat: pl.DataFrame) -> pl.DataFrame:
    """Aggregate subset of profiles for each compound."""
    feat_cols = [i for i in dat.columns if "Metadata" not in i]

    if method == "all":
        agg_df = pl.from_pandas(
            pycytominer.aggregate(
                dat.to_pandas(),
                strata=["Metadata_OASIS_ID"],
                features=feat_cols,
            ),
        )

    elif method == "allpod":
        agg_df = pl.from_pandas(
            pycytominer.aggregate(
                dat.filter(
                    pl.col("Metadata_Log10Conc") > pl.col("Metadata_POD"),
                ).to_pandas(),
                strata=["Metadata_OASIS_ID"],
                features=feat_cols,
            ),
        )

        # if no pod, then use all profiles
        dat_filt = dat.filter(
                ~pl.col("Metadata_OASIS_ID").is_in(agg_df.select(pl.col("Metadata_OASIS_ID")).to_series().to_list())
            )
        agg_df = pl.concat([
            agg_df,
            pl.from_pandas(
            pycytominer.aggregate(
                    dat_filt.to_pandas(),
                    strata=["Metadata_OASIS_ID"],
                    features=feat_cols,
                ),
            )
        ], how="vertical")

    elif method == "allpodcc":
        agg_df = pl.from_pandas(
            pycytominer.aggregate(
                dat.filter(
                    (pl.col("Metadata_Log10Conc") > pl.col("Metadata_POD"))
                    & (pl.col("Metadata_Log10Conc") < pl.col("Metadata_ccPOD")),
                ).to_pandas(),
                strata=["Metadata_OASIS_ID"],
                features=feat_cols,
            ),
        )

        # if no allpodcc, then use first profile after pod
        dat_filt = dat.filter(
                ~pl.col("Metadata_OASIS_ID").is_in(agg_df.select(pl.col("Metadata_OASIS_ID")).to_series().to_list())
            )
        agg_df = pl.concat([
            agg_df,
            pl.from_pandas(
            pycytominer.aggregate(
                    dat_filt.filter(
                        pl.col("Metadata_Concentration") == pl.col("Metadata_MinConc"),
                    ).to_pandas(),
                    strata=["Metadata_OASIS_ID"],
                    features=feat_cols,
                ),
            )
        ], how="vertical")

        # if still no pod, then use all
        dat_filt = dat.filter(
                ~pl.col("Metadata_OASIS_ID").is_in(agg_df.select(pl.col("Metadata_OASIS_ID")).to_series().to_list())
            )
        agg_df = pl.concat([
            agg_df,
            pl.from_pandas(
            pycytominer.aggregate(
                    dat_filt.to_pandas(),
                    strata=["Metadata_OASIS_ID"],
                    features=feat_cols,
                ),
            )
        ], how="vertical")

    # Annotate with aggregation and feature type
    agg_df = agg_df.with_columns(pl.lit(method).alias("Metadata_AggType"))

    return agg_df


def aggregate_profiles(
    prof_path: str,
    pod_path: str,
    agg_path: str,
) -> None:
    """Aggregate subset of profiles for each compound."""
    # 1. Read in data
    profiles = pl.read_parquet(prof_path).rename({"Metadata_Count_Cells": "Cell_Count"})
    pods = pl.read_parquet(pod_path)

    # 2. Process metadata

    # remove controls
    controls = ["DMSO"]
    profiles = profiles.filter(~pl.col("Metadata_Compound").is_in(controls))

    # Add POD
    profiles = profiles.join(
        pods.select(["Metadata_Compound", "bmd", "cc_POD"]).rename({
            "bmd": "Metadata_POD",
            "cc_POD": "Metadata_ccPOD",
        }),
        on="Metadata_Compound",
        how="left",
    )

    # Identify first conc after POD and last conc below ccPOD
    min_conc = (
        profiles.filter(pl.col("Metadata_Log10Conc") > pl.col("Metadata_POD"))
        .group_by(["Metadata_OASIS_ID"])
        .agg(pl.min("Metadata_Concentration").alias("Metadata_MinConc"))
    )

    max_conc = (
        profiles.filter(
            (pl.col("Metadata_Log10Conc") > pl.col("Metadata_POD"))
            & (pl.col("Metadata_Log10Conc") < pl.col("Metadata_ccPOD")),
        )
        .group_by(["Metadata_OASIS_ID"])
        .agg(pl.max("Metadata_Concentration").alias("Metadata_MaxConc"))
    )

    profiles = profiles.join(min_conc, on="Metadata_OASIS_ID", how="left").join(
        max_conc,
        on="Metadata_OASIS_ID",
        how="left",
    )

    # 3. Aggregate profiles
    methods = ["all", "allpod", "allpodcc"]
    agg_df = []

    for method in methods:
        agg_df.append(aggregate_compound(method, profiles))

    agg_df = pl.concat(agg_df)

    # 4. Write out results
    agg_df.write_parquet(agg_path)
