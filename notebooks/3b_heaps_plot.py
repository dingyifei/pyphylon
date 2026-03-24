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

    from pyphylon.pangenome import estimate_pan_core_size, fit_heaps_by_iteration

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    mo.md(r"""
    # 3b: Heaps' Law — Pangenome Growth Analysis

    Heaps' law models pangenome growth as genomes are accumulated:

    $$G(n) = \kappa \cdot n^\gamma$$

    where $\kappa$ is a scaling constant and $\gamma$ (gamma/lambda) controls
    the growth rate. A **gamma < 1** indicates an **open pangenome** — new
    genes continue to appear as more genomes are added, suggesting the
    species has access to a large and diverse gene pool.

    This notebook estimates pan, core, accessory, and rare genome size
    curves by random genome accumulation, then fits Heaps' law to each
    segment to produce a stacked growth plot.
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Setup
    """)
    return


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
    return DATA, FIG, OUT, SPECIES, TEMP


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
    return (df_genes,)


@app.cell
def _(df_genes):
    """Estimate pan/core/acc/rare genome size curves (cached — ~20 min)."""
    df_genes_sparse = df_genes.astype(pd.SparseDtype("int8", 0))

    with mo.persistent_cache("3b_pan_core"):
        df_pan_core = estimate_pan_core_size(df_genes_sparse, num_iter=20, log_batch=1)

    mo.output.replace(
        mo.md(
            f"Pan/core estimation complete: **{df_pan_core.shape[0]}** iterations x **{df_pan_core.shape[1]}** columns"
        )
    )
    return (df_pan_core,)


@app.cell
def _(OUT, df_genes, df_pan_core):
    """Fit Heaps' law to each CAR section and save summary."""
    output_pan = fit_heaps_by_iteration(df_pan_core, section="pan")
    output_core = fit_heaps_by_iteration(df_pan_core, section="core")
    output_acc = fit_heaps_by_iteration(df_pan_core, section="acc")
    output_rare = fit_heaps_by_iteration(df_pan_core, section="rare")

    # NOTE: If this errors with KeyError (e.g. 'Pan465'), the persistent cache
    # "3b_pan_core" is stale — delete it and re-run so df_pan_core matches df_genes.
    n_genomes = df_genes.shape[1]
    last_col = str(n_genomes)

    summary = pd.DataFrame(
        {
            "section": ["pan", "core", "accessory", "rare"],
            "lambda_mean": [
                output_pan.lambda_.mean(),
                output_core.lambda_.mean(),
                output_acc.lambda_.mean(),
                output_rare.lambda_.mean(),
            ],
            "kappa_mean": [
                output_pan.kappa.mean(),
                output_core.kappa.mean(),
                output_acc.kappa.mean(),
                output_rare.kappa.mean(),
            ],
            "gene_count_at_n": [
                df_pan_core[f"Pan{last_col}"].mean(),
                df_pan_core[f"Core{last_col}"].mean(),
                df_pan_core[f"Acc{last_col}"].mean(),
                df_pan_core[f"Rare{last_col}"].mean(),
            ],
            "n_genomes": [n_genomes] * 4,
            "n_iterations": [df_pan_core.shape[0]] * 4,
        }
    )
    summary.to_csv(os.path.join(OUT, "data", "3b_heaps_summary.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"Heaps' law fit ({n_genomes} genomes, {df_pan_core.shape[0]} iterations):\n\n"
            f"| Section | Lambda | Kappa | Genes at N |\n"
            f"|---------|--------|-------|------------|\n"
            f"| Pan | {output_pan.lambda_.mean():.4f} | {output_pan.kappa.mean():.1f} | {df_pan_core[f'Pan{last_col}'].mean():,.0f} |\n"
            f"| Core | {output_core.lambda_.mean():.4f} | {output_core.kappa.mean():.1f} | {df_pan_core[f'Core{last_col}'].mean():,.0f} |\n"
            f"| Accessory | {output_acc.lambda_.mean():.4f} | {output_acc.kappa.mean():.1f} | {df_pan_core[f'Acc{last_col}'].mean():,.0f} |\n"
            f"| Rare | {output_rare.lambda_.mean():.4f} | {output_rare.kappa.mean():.1f} | {df_pan_core[f'Rare{last_col}'].mean():,.0f} |"
        )
    )
    return output_acc, output_core, output_rare


@app.cell
def _():
    mo.md("""
    ## Heaps' Law Stacked Plot

    The stacked area plot below shows fitted growth curves for the core,
    accessory, and rare genome segments on a log-linear scale. The total
    height at any point represents the predicted pangenome size for that
    number of genomes.
    """)
    return


@app.cell
def _(FIG, df_genes, output_acc, output_core, output_rare):
    """Create and save Heaps' law stacked area plot."""
    x = list(range(1, df_genes.shape[1] + 1))
    y_core = output_core.kappa.mean() * np.array(x) ** output_core.lambda_.mean()
    y_acc = output_acc.kappa.mean() * np.array(x) ** output_acc.lambda_.mean()
    y_rare = output_rare.kappa.mean() * np.array(x) ** output_rare.lambda_.mean()

    fig_heaps, ax_heaps = plt.subplots(figsize=(8, 5))
    ax_heaps.stackplot(x, y_core, y_acc, y_rare, labels=["Core", "Accessory", "Rare"])
    ax_heaps.set_yscale("log")
    ax_heaps.set_xlabel("Number of genomes")
    ax_heaps.set_ylabel("Number of gene clusters")
    ax_heaps.set_title("Heaps' Law (Log-linear)")
    ax_heaps.legend(loc="upper left")
    ax_heaps.grid(False)
    plt.tight_layout()

    fig_heaps.savefig(os.path.join(FIG, "3b_heaps_law.png"), bbox_inches="tight")
    mo.output.replace(fig_heaps)
    return


if __name__ == "__main__":
    app.run()
