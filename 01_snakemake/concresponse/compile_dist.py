import polars as pl


def filter_dist(thresh: int, dist: pl.DataFrame) -> pl.DataFrame:
    meta_cols = [i for i in dist.columns if "Metadata" in i]
    dist_meta = dist.select(meta_cols)

    dist = dist.with_columns(
        pl.when(pl.col("Metadata_well_type") == "DMSO")
        .then(
            pl.concat_str(
                [pl.col("Metadata_Plate"), pl.col("Metadata_Compound")], separator="_"
            )
        )
        .otherwise(
            pl.concat_str(
                [pl.col("Metadata_Compound"), pl.col("Metadata_Concentration")],
                separator="_",
            )
        )
        .alias("Metadata_Group")
    )

    # Get column types
    dist_cols = [i for i in dist.columns if "Metadata" not in i]

    dist_long = dist.unpivot(
        index=["Metadata_Group", "Metadata_Plate", "Metadata_Well"],
        on=dist_cols,
        variable_name="Distance_type",
        value_name="Distance",
    )

    # Compute median and mad for each sample group
    dist_median = (
        dist.group_by("Metadata_Group")
        .agg([pl.median(col).alias(col) for col in dist_cols])
        .unpivot(
            index="Metadata_Group",
            on=dist_cols,
            variable_name="Distance_type",
            value_name="Median_distance",
        )
    )

    dist_mad = (
        dist.group_by("Metadata_Group")
        .agg(
            [
                pl.col(col)
                .map_elements(
                    lambda x: (x - x.median()).abs().median(), return_dtype=pl.Float64
                )
                .alias(col)
                for col in dist_cols
            ]
        )
        .unpivot(
            index="Metadata_Group",
            on=dist_cols,
            variable_name="Distance_type",
            value_name="MAD",
        )
    )

    # Compute deviation relative to MAD and filter
    dist_long = (
        dist_long.join(dist_median, on=["Metadata_Group", "Distance_type"])
        .join(dist_mad, on=["Metadata_Group", "Distance_type"])
        .with_columns(
            (pl.col("Distance") - pl.col("Median_distance"))
            .abs()
            .alias("Absolute_Deviation")
        )
        .with_columns(
            (pl.col("Absolute_Deviation") / pl.col("MAD"))
            .log(base=2)
            .alias("Log2_AD_MAD")
        )
        .with_columns(
            pl.when(pl.col("Absolute_Deviation") == 0)
            .then(pl.col("Distance"))
            .when(pl.col("Log2_AD_MAD").abs() > thresh)
            .then(pl.lit(None))
            .otherwise(pl.col("Distance"))
            .alias("Distance_filtered")
        )
    )

    # Convert back to wide format
    dist_filt = (
        dist_long.select(
            ["Metadata_Plate", "Metadata_Well", "Distance_type", "Distance_filtered"]
        )
        .pivot(
            index=["Metadata_Plate", "Metadata_Well"],
            on="Distance_type",
            values="Distance_filtered",
        )
        .drop_nulls()
    )

    # put metadata back
    dist_filt = dist_meta.join(dist_filt, on=["Metadata_Plate", "Metadata_Well"])

    return dist_filt


def compile_dist(
    input_files: list, transform: str, thresh: int, output_path: str
) -> None:
    dfs = []
    for fp in input_files:
        dat = pl.read_parquet(fp)
        dfs.append(dat.unique())

    dfs = pl.concat(dfs, how="vertical")
    meta_cols = [i for i in dfs.columns if "Metadata" in i and i != "Metadata_Distance"]
    df_wide = dfs.pivot(
        values="Distance",
        index=meta_cols,
        columns="Metadata_Distance",
        aggregate_function="median",
    )

    # apply transform if specified
    if transform == "log10":
        dist_cols = [i for i in df_wide.columns if "Metadata" not in i]
        df_wide = df_wide.with_columns(
            [pl.col(col).log10().alias(col) for col in dist_cols]
        )

        # shift by min_val to ensure all log10 transformed values are > 0
        min_val = (
            df_wide.select([pl.col(col).min() for col in dist_cols])
            .transpose()
            .min()
            .item()
        )
        df_wide = df_wide.with_columns(
            [(pl.col(col) - min_val + 1e-6).alias(col) for col in dist_cols]
        )

    # Filter out outlier replicates
    df_wide_filt = filter_dist(thresh, df_wide)

    df_wide_filt.write_parquet(output_path)
