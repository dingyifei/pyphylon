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
    from sklearn.cluster import KMeans
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score

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

    TEMP = CONFIG["TEMP_DIR"]
    DATA = CONFIG["DATA_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    os.makedirs(os.path.join(DATA, "processed", "nmf-outputs"), exist_ok=True)

    return CONFIG, DATA, FIG, OUT, TEMP


@app.cell
def _(DATA, TEMP):
    """Load normalized NMF matrices and accessory genome data."""
    # L and A are already normalized by 4a (99th percentile normalization)
    L_norm = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "L.csv"), index_col=0)
    A_norm = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "A.csv"), index_col=0)

    # Rename numeric indices to phylonN
    L_norm.columns = [f"phylon{c}" for c in L_norm.columns]
    A_norm.index = [f"phylon{i}" for i in A_norm.index]

    # Accessory genome for gene frequency calculation
    df_acc = pd.read_csv(os.path.join(DATA, "processed", "CAR_genomes", "df_acc.csv"), index_col=0)

    # Metadata (for display/reporting)
    metadata = pd.read_csv(os.path.join(TEMP, "2d_enriched_metadata.csv"), dtype={"genome_id": str})

    n_phylons = L_norm.shape[1]
    mo.output.replace(
        mo.md(
            f"Loaded inputs:\n\n"
            f"- **L_norm:** {L_norm.shape[0]} genes x {L_norm.shape[1]} phylons (already normalized)\n"
            f"- **A_norm:** {A_norm.shape[0]} phylons x {A_norm.shape[1]} genomes (already normalized)\n"
            f"- **df_acc:** {df_acc.shape[0]} genes x {df_acc.shape[1]} genomes\n"
            f"- **Metadata:** {metadata.shape[0]} genomes\n"
            f"- **Number of phylons:** {n_phylons}"
        )
    )
    return L_norm, A_norm, df_acc, metadata, n_phylons


@app.cell
def _(L_norm):
    """Binarize L_norm via 3-means clustering (highest-center cluster per column)."""
    L_binarized = pd.DataFrame(0, index=L_norm.index, columns=L_norm.columns)

    for col in L_norm.columns:
        km = KMeans(n_clusters=3, random_state=0, n_init="auto")
        km.fit(L_norm[col].values.reshape(-1, 1))
        highest = np.argmax(km.cluster_centers_)
        L_binarized[col] = (km.labels_ == highest).astype(int)

    genes_per_phylon = L_binarized.sum(axis=0)
    mo.output.replace(
        mo.md(
            f"L binarized: {L_binarized.shape[0]} genes x {L_binarized.shape[1]} phylons\n\n"
            f"Genes per phylon:\n\n"
            + "\n".join(f"- **{col}:** {genes_per_phylon[col]} genes" for col in L_binarized.columns)
        )
    )
    return (L_binarized,)


@app.cell
def _(A_norm):
    """Binarize A_norm via fixed threshold (>= 0.5)."""
    A_binarized = (A_norm >= 0.5).astype(int)

    strains_per_phylon = A_binarized.sum(axis=1)
    mo.output.replace(
        mo.md(
            f"A binarized: {A_binarized.shape[0]} phylons x {A_binarized.shape[1]} genomes\n\n"
            f"Strains per phylon:\n\n"
            + "\n".join(f"- **{idx}:** {strains_per_phylon[idx]} strains" for idx in A_binarized.index)
        )
    )
    return (A_binarized,)


@app.cell
def _(A_binarized, L_binarized):
    """Build sorted gene, phylon, and strain orderings for visualization."""
    # --- Phylon order from L_binarized dendrogram ---
    g_bin = sns.clustermap(
        L_binarized,
        method="ward",
        metric="euclidean",
        cmap="Greys",
        yticklabels=False,
    )
    phylon_order = g_bin.data2d.columns.tolist()
    plt.close(g_bin.fig)

    # --- Gene order: zero-phylon -> single-phylon (per phylon) -> poly-phylon (subclustered) ---
    gene_order = []

    # Zero-phylon genes
    zero_cond = L_binarized.sum(axis=1) == 0
    gene_order.extend(L_binarized[zero_cond].index.tolist())

    # Single-phylon genes (ordered by phylon_order)
    for phylon in phylon_order:
        single_cond = L_binarized.sum(axis=1) == 1
        in_phylon = L_binarized[phylon] == 1
        gene_order.extend(L_binarized[in_phylon & single_cond].index.tolist())

    # Poly-phylon genes (subclustered by ward linkage)
    max_n_phylons = int(L_binarized.sum(axis=1).max())
    for n_active in range(2, max_n_phylons + 1):
        num_cond = L_binarized.sum(axis=1) == n_active
        subset = L_binarized[num_cond]
        if subset.empty:
            continue
        if len(subset) == 1:
            gene_order.extend(subset.index.tolist())
            continue
        gg = sns.clustermap(subset, method="ward", metric="euclidean", col_cluster=False, yticklabels=False)
        gene_order.extend(gg.data2d.index.tolist())
        plt.close(gg.fig)

    # --- Strain order: no-phylon -> single-phylon (per phylon) -> multi-phylon ---
    strain_order = []
    assigned_strains = set()

    # No-phylon strains
    no_phylon_mask = A_binarized.sum(axis=0) == 0
    no_phylon_strains = A_binarized.columns[no_phylon_mask].tolist()
    strain_order.extend(no_phylon_strains)
    assigned_strains.update(no_phylon_strains)

    single_phylon_strains = A_binarized.columns[A_binarized.sum(axis=0) == 1]
    multi_phylon_strains = A_binarized.columns[A_binarized.sum(axis=0) > 1]

    # Single then multi, per phylon
    for phylon in phylon_order:
        # Single-phylon strains affiliated to this phylon
        in_single = A_binarized.loc[phylon, single_phylon_strains] == 1
        new_single = [s for s in single_phylon_strains[in_single] if s not in assigned_strains]
        strain_order.extend(new_single)
        assigned_strains.update(new_single)

        # Multi-phylon strains affiliated to this phylon
        in_multi = A_binarized.loc[phylon, multi_phylon_strains] == 1
        new_multi = [s for s in multi_phylon_strains[in_multi] if s not in assigned_strains]
        strain_order.extend(new_multi)
        assigned_strains.update(new_multi)

    # Any remaining strains not yet assigned
    remaining = [s for s in A_binarized.columns if s not in assigned_strains]
    strain_order.extend(remaining)

    mo.output.replace(
        mo.md(
            f"Orderings built:\n\n"
            f"- **Phylon order:** {phylon_order}\n"
            f"- **Gene order:** {len(gene_order)} genes\n"
            f"- **Strain order:** {len(strain_order)} strains"
        )
    )
    return gene_order, phylon_order, strain_order


@app.cell
def _(FIG, L_binarized, df_acc):
    """Regression: gene frequency vs number of active phylons per gene."""
    gene_freq = df_acc.sum(axis=1).reindex(L_binarized.index)
    n_phylons_per_gene = L_binarized.sum(axis=1)

    model = LinearRegression()
    X_reg = gene_freq.values.reshape(-1, 1)
    y_reg = n_phylons_per_gene.values
    model.fit(X_reg, y_reg)
    y_pred = model.predict(X_reg)

    slope = model.coef_[0]
    intercept = model.intercept_
    r2 = r2_score(y_reg, y_pred)

    fig_reg, ax_reg = plt.subplots(figsize=(8, 5))
    sns.regplot(x=gene_freq, y=n_phylons_per_gene, ax=ax_reg)
    ax_reg.set_xlabel("Gene frequency (number of genomes)")
    ax_reg.set_ylabel("Number of active phylons")
    ax_reg.set_title(f"Gene Frequency vs Phylon Count (R² = {r2:.4f})")
    plt.tight_layout()

    fig_reg.savefig(os.path.join(FIG, "5a_regression.png"), bbox_inches="tight")
    mo.output.replace(
        mo.vstack(
            [
                mo.md(f"Linear regression: y = {slope:.6f}x + {intercept:.4f}\n\nR² = {r2:.4f}"),
                fig_reg,
            ]
        )
    )
    return (r2,)


@app.cell
def _(A_binarized, FIG, L_binarized, gene_order, phylon_order, strain_order):
    """Generate sorted heatmaps for L_binarized and A_binarized."""
    # L_binarized: genes (rows, sorted) x phylons (columns, sorted)
    g_L = sns.clustermap(
        L_binarized.loc[gene_order, phylon_order],
        row_cluster=False,
        col_cluster=False,
        yticklabels=False,
        cmap="Greys",
        figsize=(6, 12),
    )
    g_L.fig.suptitle("L_binarized (sorted genes x phylons)", y=1.02)
    g_L.savefig(os.path.join(FIG, "5a_L_binarized_sorted.png"), bbox_inches="tight")

    # A_binarized: phylons (rows, sorted) x strains (columns, sorted)
    g_A = sns.clustermap(
        A_binarized.loc[phylon_order, strain_order],
        row_cluster=False,
        col_cluster=False,
        xticklabels=False,
        cmap="Greys",
        figsize=(12, 6),
    )
    g_A.fig.suptitle("A_binarized (sorted phylons x strains)", y=1.02)
    g_A.savefig(os.path.join(FIG, "5a_A_binarized_sorted.png"), bbox_inches="tight")

    mo.output.replace(mo.vstack([g_L.fig, g_A.fig]))


@app.cell
def _(A_binarized, DATA, L_binarized, L_norm, A_norm, OUT, n_phylons, r2):
    """Save normalized and binarized matrices, plus summary CSV."""
    nmf_dir = os.path.join(DATA, "processed", "nmf-outputs")

    # Save normalized matrices (copies of L.csv/A.csv, for downstream compatibility)
    L_norm.to_csv(os.path.join(nmf_dir, "L_norm.csv"))
    A_norm.to_csv(os.path.join(nmf_dir, "A_norm.csv"))

    # Save binarized matrices (overwrite 4a's versions — 5a is canonical owner)
    L_binarized.to_csv(os.path.join(nmf_dir, "L_binarized.csv"))
    A_binarized.to_csv(os.path.join(nmf_dir, "A_binarized.csv"))

    # Summary CSV for Quarto report
    summary = pd.DataFrame(
        {
            "metric": [
                "n_phylons",
                "n_genes",
                "n_genomes",
                "regression_r2",
                "zero_phylon_genes",
                "single_phylon_genes",
                "multi_phylon_genes",
            ],
            "value": [
                n_phylons,
                L_binarized.shape[0],
                A_binarized.shape[1],
                round(r2, 4),
                int((L_binarized.sum(axis=1) == 0).sum()),
                int((L_binarized.sum(axis=1) == 1).sum()),
                int((L_binarized.sum(axis=1) > 1).sum()),
            ],
        }
    )
    summary.to_csv(os.path.join(OUT, "data", "5a_phylon_summary.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"Saved outputs:\n\n"
            f"- `{nmf_dir}/L_norm.csv` ({L_norm.shape})\n"
            f"- `{nmf_dir}/A_norm.csv` ({A_norm.shape})\n"
            f"- `{nmf_dir}/L_binarized.csv` ({L_binarized.shape})\n"
            f"- `{nmf_dir}/A_binarized.csv` ({A_binarized.shape})\n"
            f"- Summary: {n_phylons} phylons, R² = {r2:.4f}"
        )
    )


if __name__ == "__main__":
    app.run()
