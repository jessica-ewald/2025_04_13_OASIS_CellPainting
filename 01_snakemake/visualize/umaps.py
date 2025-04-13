import anndata
import matplotlib
import matplotlib.pyplot as plt
import polars as pl
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use("Agg")


def make_umaps(prof_path: str, morph_pod: str, cc_pod: str, ldh_pod: str, mtt_pod: str, plot_path: str) -> None:
    data = pl.read_parquet(prof_path)

    cc = (
        pl.read_parquet(cc_pod)
        .filter(
            pl.col("all.pass") == True,
        )
        .select(["Metadata_Compound", "bmd"])
        .rename({"bmd": "Metadata_cc_POD"})
    )
    ldh = (
        pl.read_parquet(ldh_pod)
        .filter(
            pl.col("all.pass") == True,
        )
        .select(["Metadata_Compound", "bmd"])
        .rename({"bmd": "Metadata_ldh_POD"})
    )
    mtt = (
        pl.read_parquet(mtt_pod)
        .filter(
            pl.col("all.pass") == True,
        )
        .select(["Metadata_Compound", "bmd"])
        .rename({"bmd": "Metadata_mtt_POD"})
    )
    morph = pl.read_parquet(morph_pod).select(["Metadata_Compound", "bmd"]).rename({"bmd": "Metadata_morph_POD"})

    # Add PODs to metadata
    data = data.join(cc, on="Metadata_Compound", how="left")
    data = data.join(ldh, on="Metadata_Compound", how="left")
    data = data.join(mtt, on="Metadata_Compound", how="left")
    data = data.join(morph, on="Metadata_Compound", how="left")

    # Add columns to label different sample subsets based on their bioactivity
    data = data.with_columns(
        (pl.col("Metadata_Log10Conc") > pl.col("Metadata_morph_POD")).alias("Metadata_Bioactive"),
        (pl.col("Metadata_Log10Conc") < pl.col("Metadata_cc_POD")).alias("Metadata_No_Cytotox"),
    )

    data = data.with_columns(
        pl.when(pl.col("Metadata_Bioactive") == False)
        .then(False)
        .when(pl.col("Metadata_Bioactive") == True)
        .then(True)
        .otherwise(False)
        .alias("Metadata_Bioactive"),
    )
    data = data.with_columns(
        pl.when(pl.col("Metadata_No_Cytotox") == False)
        .then(False)
        .when(pl.col("Metadata_No_Cytotox") == True)
        .then(True)
        .otherwise(True)
        .alias("Metadata_No_Cytotox"),
    )

    metadata_cols = [col for col in data.columns if "Metadata" in col]

    data = data.to_pandas()
    data.sort_values(["Metadata_Plate", "Metadata_Well"], inplace=True)
    data.index = [f"{row['Metadata_Plate']}__{row['Metadata_Well']}" for _, row in data.iterrows()]
    data = data.loc[~data.index.duplicated(keep="first")]

    metadata = data[metadata_cols]
    adata = anndata.AnnData(X=data.drop(metadata_cols, axis=1))
    adata.obs = adata.obs.merge(metadata, left_index=True, right_index=True)

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Plot UMAPs (only bioactive samples)
    bioactive_mask = adata.obs["Metadata_Bioactive"]
    adata_bioactive = adata[bioactive_mask].copy()

    sc.pp.neighbors(adata_bioactive)
    sc.tl.umap(adata_bioactive)

    # Plot UMAPs (only bioactive & non-cytotoxic samples)
    nocytotox_mask = adata_bioactive.obs["Metadata_No_Cytotox"]
    adata_nocytotox = adata_bioactive[nocytotox_mask].copy()

    sc.pp.neighbors(adata_nocytotox)
    sc.tl.umap(adata_nocytotox)

    with PdfPages(plot_path) as pdf:
        sc.pl.embedding(
            adata,
            "X_umap",
            color="Metadata_source",
            s=10,
            show=False,
            title="All samples (source)",
        )
        pdf.savefig()
        plt.close()

        sc.pl.embedding(
            adata,
            "X_umap",
            color="Metadata_Count_Cells",
            s=10,
            show=False,
            title="All samples (cell count)",
        )
        pdf.savefig()
        plt.close()

        sc.pl.embedding(
            adata_bioactive,
            "X_umap",
            color="Metadata_Count_Cells",
            s=10,
            show=False,
            title="Bioactive samples (cell count)",
        )
        pdf.savefig()
        plt.close()

        sc.pl.embedding(
            adata_nocytotox,
            "X_umap",
            color="Metadata_Count_Cells",
            s=10,
            show=False,
            title="Bioactive and non-cytotoxic samples (cell count)",
        )
        pdf.savefig()
        plt.close()
