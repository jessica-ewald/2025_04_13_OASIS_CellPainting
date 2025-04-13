import polars as pl


def compile_dist(input_files: list, transform: str, output_path: str) -> None:
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
        df_wide = df_wide.with_columns([pl.col(col).log10().alias(col) for col in dist_cols])

        # shift by min_val to ensure all log10 transformed values are > 0
        min_val = df_wide.select([pl.col(col).min() for col in dist_cols]).transpose().min().item()
        df_wide = df_wide.with_columns([(pl.col(col) - min_val + 1e-6).alias(col) for col in dist_cols])

    df_wide.write_parquet(output_path)
