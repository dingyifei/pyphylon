import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml

    from pyphylon.pangenome import find_pangenome_segments

    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    """Parse config and set up directories."""
    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    TEMP = CONFIG["TEMP_DIR"]
    DATA = CONFIG["DATA_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]
    SPECIES = CONFIG["PG_NAME"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    os.makedirs(os.path.join(DATA, "processed", "CAR_genomes"), exist_ok=True)

    return CONFIG, DATA, FIG, OUT, SPECIES, TEMP


@app.cell
def _(DATA, SPECIES, TEMP):
    """Load P matrix (gene presence/absence) and enriched metadata."""
    df_genes = pd.read_pickle(os.path.join(DATA, "processed", "cd-hit-results", f"{SPECIES}_strain_by_gene.pickle.gz"))
    df_genes = df_genes.fillna(0).sparse.to_dense().astype("int8")

    metadata = pd.read_csv(
        os.path.join(TEMP, "2d_enriched_metadata.csv"),
        dtype={"genome_id": str},
    )
    metadata = metadata[metadata["genome_id"].isin(df_genes.columns)]

    mo.output.replace(
        mo.md(
            f"Loaded inputs:\n\n"
            f"- **P matrix:** {df_genes.shape[0]} genes x {df_genes.shape[1]} genomes\n"
            f"- **Metadata:** {metadata.shape[0]} genomes (filtered to P matrix)"
        )
    )
    return df_genes, metadata


@app.cell
def _(FIG, df_genes):
    """Plot gene frequency distribution across all genomes."""
    df_gene_freq = df_genes.sum(axis=1)

    fig_freq, ax_freq = plt.subplots(figsize=(8, 4))
    sns.histplot(df_gene_freq, binwidth=50, ax=ax_freq)
    ax_freq.set_yscale("log")
    ax_freq.set_xlabel("Gene frequency (number of genomes)")
    ax_freq.set_ylabel("Number of genes")
    ax_freq.set_title("Gene Frequency Distribution")
    plt.tight_layout()

    fig_freq.savefig(os.path.join(FIG, "3a_gene_frequency.png"), bbox_inches="tight")
    mo.output.replace(fig_freq)


@app.cell
def _(FIG, df_genes):
    """Fit pangenome frequency model and identify core/rare boundaries."""
    fig_seg, ax_seg = plt.subplots()
    segments, _popt, r_squared, mae, _ax = find_pangenome_segments(df_genes, threshold=0.1, ax=ax_seg)
    fig_seg.savefig(os.path.join(FIG, "3a_pangenome_segments.png"), bbox_inches="tight")

    mo.output.replace(
        mo.md(
            f"Pangenome segment analysis:\n\n"
            f"- **Core threshold:** freq > {np.floor(segments[0]):.0f} genomes\n"
            f"- **Rare threshold:** freq < {np.ceil(segments[1]):.0f} genomes\n"
            f"- **R²:** {r_squared:.4f}\n"
            f"- **MAE:** {mae:.4f}"
        )
    )
    return (segments,)


@app.cell
def _(DATA, OUT, df_genes, segments):
    """Classify genes into Core/Accessory/Rare and save results."""
    df_freq = df_genes.sum(axis=1)

    df_core = df_genes[df_freq > np.floor(segments[0])]
    df_rare = df_genes[df_freq < np.ceil(segments[1])]
    acc_gene_list = list(set(df_genes.index) - set(df_core.index) - set(df_rare.index))
    df_acc = df_genes.loc[acc_gene_list].copy()

    # Save CAR matrices (index=True: gene IDs are the row index, downstream reads with index_col=0)
    car_dir = os.path.join(DATA, "processed", "CAR_genomes")
    df_core.to_csv(os.path.join(car_dir, "df_core.csv"))
    df_acc.to_csv(os.path.join(car_dir, "df_acc.csv"))
    df_rare.to_csv(os.path.join(car_dir, "df_rare.csv"))

    # Save summary for Quarto report
    _summary = pd.DataFrame(
        {
            "category": ["Core", "Accessory", "Rare", "Total"],
            "n_genes": [df_core.shape[0], df_acc.shape[0], df_rare.shape[0], df_genes.shape[0]],
            "n_genomes": [df_core.shape[1], df_acc.shape[1], df_rare.shape[1], df_genes.shape[1]],
            "threshold": [
                f"> {np.floor(segments[0]):.0f}",
                f"{np.ceil(segments[1]):.0f} - {np.floor(segments[0]):.0f}",
                f"< {np.ceil(segments[1]):.0f}",
                "all",
            ],
        }
    )
    _summary.to_csv(os.path.join(OUT, "data", "3a_car_summary.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"CAR classification saved:\n\n"
            f"| Category | Genes | Threshold |\n"
            f"|----------|-------|-----------|\n"
            f"| Core | {df_core.shape[0]} | > {np.floor(segments[0]):.0f} |\n"
            f"| Accessory | {df_acc.shape[0]} | {np.ceil(segments[1]):.0f} - {np.floor(segments[0]):.0f} |\n"
            f"| Rare | {df_rare.shape[0]} | < {np.ceil(segments[1]):.0f} |\n"
            f"| **Total** | **{df_genes.shape[0]}** | |"
        )
    )


if __name__ == "__main__":
    app.run()
