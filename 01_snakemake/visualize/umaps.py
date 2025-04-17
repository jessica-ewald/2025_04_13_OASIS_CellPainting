import anndata
import matplotlib
import matplotlib.pyplot as plt
import polars as pl
import scanpy as sc
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use("Agg")


def make_umaps(prof_path: str, morph_pod: str, cc_pod: str, plot_path: str) -> None:
    data = pl.read_parquet(prof_path)

    cc = (
        pl.read_parquet(cc_pod)
        .filter((pl.col("all.pass") == True))
        .select(["Metadata_Compound", "bmd"])
        .rename({"bmd": "Metadata_cc_POD"})
    )
    morph = (
        pl.read_parquet(morph_pod)
        .select(["Metadata_Compound", "bmd"])
        .rename({"bmd": "Metadata_morph_POD"})
    )

    # Add PODs to metadata
    data = data.join(cc, on="Metadata_Compound", how="left")
    data = data.join(morph, on="Metadata_Compound", how="left")

    # Add columns to label different sample subsets based on their bioactivity
    data = data.with_columns(
        (pl.col("Metadata_Log10Conc") > pl.col("Metadata_morph_POD")).alias(
            "Metadata_Bioactive"
        ),
        (pl.col("Metadata_Log10Conc") < pl.col("Metadata_cc_POD")).alias(
            "Metadata_No_Cytotox"
        ),
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
    ).sample(fraction=1.0, seed=42, shuffle=True)

    metadata_cols = [col for col in data.columns if "Metadata" in col]

    data = data.to_pandas()
    data.index = [
        f"{row['Metadata_Plate']}__{row['Metadata_Well']}" for _, row in data.iterrows()
    ]
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
            color="Metadata_Count_Cells",
            s=10,
            show=False,
            title="All samples (cell count)",
        )
        pdf.savefig()
        plt.close()

        sc.pl.embedding(
            adata,
            "X_umap",
            color="Metadata_Source",
            s=10,
            show=False,
            title="All samples (batch)",
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

        # Add cell count boxplot
        plt.figure(figsize=(10, 6))
        sns.boxplot(
            data=data[data["Metadata_well_type"] == "DMSO"],
            x="Metadata_Source",
            y="Metadata_Count_Cells",
            palette="pastel"
        )
        plt.xlabel("Metadata Source")
        plt.ylabel("Cell Count")
        plt.title("DMSO Cell Count by Source")
        plt.tight_layout()

        pdf.savefig()
        plt.close()

        # Add cell count boxplot by plate
        sns.set(style="whitegrid")
        g = sns.catplot(
            data=data[data["Metadata_well_type"] == "DMSO"],
            x="Metadata_Plate",
            y="Metadata_Count_Cells",
            col="Metadata_Source",
            kind="box",
            col_wrap=2,
            height=4,
            aspect=1.5,
            palette="pastel",
            sharex=False,
            sharey=True 
        )

        for ax in g.axes.flatten():
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        g.set_axis_labels("Plate", "Cell Count")
        g.fig.subplots_adjust(top=0.9)
        g.fig.suptitle("Cell Count per Plate by Batch")
        plt.tight_layout()

        pdf.savefig()
        plt.close()

        # Plot mean cell count by well position
        df = pl.read_parquet(prof_path).select(
            ["Metadata_Source", "Metadata_Plate", "Metadata_Well", "Metadata_well_type", "Metadata_Count_Cells"]
        ).group_by(["Metadata_Source", "Metadata_Well"]).agg([
            pl.col("Metadata_Count_Cells").mean()
        ]).to_pandas()

        df["Well_Row"] = df["Metadata_Well"].str.extract(r"([A-P])")
        df["Well_Col"] = df["Metadata_Well"].str.extract(r"(\d{2})")

        row_order = list("ABCDEFGHIJKLMNOP")


        sources = df["Metadata_Source"].unique()
        for source in sources:
            sub_df = df[df["Metadata_Source"] == source]
            
            heatmap_data = (
                sub_df.groupby(["Well_Row", "Well_Col"])["Metadata_Count_Cells"]
                .mean()
                .unstack()
                .reindex(index=row_order)
            )

            plt.figure(figsize=(14, 8))
            sns.heatmap(
                heatmap_data,
                cmap="viridis",
                linewidths=0.5,
                linecolor="white",
                cbar_kws={"label": "Mean Cell Count"},
            )
            plt.title(f"Cell Count by Well Position — {source}")
            plt.xlabel("Well Column")
            plt.ylabel("Well Row")
            plt.tight_layout()
            pdf.savefig()
            plt.close()

