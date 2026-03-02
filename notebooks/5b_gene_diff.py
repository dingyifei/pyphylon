import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    import yaml

    from pyphylon.plotting import generate_phylon_dendrogram

    plt.rcParams["pdf.fonttype"] = 42
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

    DATA = CONFIG["DATA_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)

    return CONFIG, DATA, FIG, OUT


@app.cell
def _(DATA):
    """Load L_binarized matrix from 5a outputs."""
    L_binarized = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "L_binarized.csv"), index_col=0)

    mo.output.replace(
        mo.md(
            f"Loaded L_binarized: **{L_binarized.shape[0]}** genes x **{L_binarized.shape[1]}** phylons\n\n"
            f"Values: {sorted(L_binarized.values.ravel().tolist()).__class__} unique = "
            f"{sorted(set(L_binarized.values.ravel().tolist()))}"
        )
    )
    return (L_binarized,)


@app.cell
def _(FIG, L_binarized):
    """Generate Ward clustermap and annotated phylon dendrogram."""
    # Ward hierarchical clustering heatmap
    g = sns.clustermap(L_binarized, method="ward", cmap="hot_r")
    g.savefig(os.path.join(FIG, "5b_clustermap.png"), bbox_inches="tight")
    plt.close(g.fig)

    # Annotated dendrogram with gene counts per split
    fig_dend, ax_dend = plt.subplots(figsize=(10, 6))
    ax_dend, df_stats, split_gene_sets_labeled = generate_phylon_dendrogram(
        L_binarized,
        text_offset=0,
        labels=["exclusive_genes", "total_split_genes"],
        orientation="left",
        ax=ax_dend,
    )
    fig_dend.tight_layout()
    fig_dend.savefig(os.path.join(FIG, "5b_phylon_dendrogram.png"), bbox_inches="tight")

    mo.output.replace(mo.vstack([mo.md("**Phylon dendrogram** (annotated with gene counts)"), fig_dend]))
    return df_stats, split_gene_sets_labeled


@app.cell
def _(OUT, df_stats, L_binarized):
    """Save dendrogram statistics."""
    df_stats.to_csv(os.path.join(OUT, "data", "5b_gene_diff_stats.csv"))

    mo.output.replace(
        mo.md(
            f"Saved outputs:\n\n"
            f"- `5b_gene_diff_stats.csv`: {df_stats.shape[0]} rows (splits + leaves)\n"
            f"- `5b_clustermap.png`: Ward clustering heatmap\n"
            f"- `5b_phylon_dendrogram.png`: annotated dendrogram\n\n"
            f"L_binarized: {L_binarized.shape[0]} genes x {L_binarized.shape[1]} phylons"
        )
    )


if __name__ == "__main__":
    app.run()
