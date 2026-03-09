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

    from pyphylon.biointerp import collect_functions, get_pg_to_locus_map
    from pyphylon.plotting_util import unique_genes_by_phylon

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    mo.md("""
    # 5a: Phylon Characterization

    Refine the NMF decomposition from step 4a into discrete phylon memberships.
    The **L matrix** (genes × phylons) is binarized via **3-means clustering** —
    for each phylon column, K-Means with k=3 is applied and genes in the
    highest-center cluster are assigned membership. The **A matrix**
    (phylons × genomes) is binarized with a fixed threshold (≥ 0.5).

    Sorted visualizations reveal the block-diagonal structure of gene–phylon
    and strain–phylon associations, and a linear regression quantifies the
    relationship between gene frequency and multi-phylon membership.
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
    os.makedirs(os.path.join(DATA, "processed", "nmf-outputs"), exist_ok=True)
    return DATA, FIG, OUT, SPECIES, TEMP


@app.cell
def _():
    mo.md("""
    ## Load Inputs
    """)
    return


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
    return A_norm, L_norm, df_acc, metadata, n_phylons


@app.cell
def _():
    mo.md("""
    ## Binarize L Matrix (3-means Clustering)

    For each phylon column in L_norm, K-Means with **k=3** is applied.
    Genes belonging to the cluster with the **highest center** are assigned
    membership (1), all others are set to 0. Using 3 clusters generally
    provides a better precision–recall tradeoff than 2-means by separating
    the high-affinity tail from intermediate and background signal.
    """)
    return


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
def _():
    mo.md("""
    ## Binarize A Matrix (Threshold)

    The A matrix (phylons × genomes) is binarized with a fixed threshold:
    strains with affinity **≥ 0.5** are assigned to the phylon. This
    simple threshold works well when the affinity distribution is bimodal.
    """)
    return


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
def _():
    mo.md("""
    ## Sorted L_binarized (Gene-Phylon Membership)

    Genes are sorted into three tiers to reveal block-diagonal structure:

    1. **Zero-phylon genes** — not assigned to any phylon
    2. **Single-phylon genes** — exclusive to one phylon (highest differentiating
       power), grouped by phylon in dendrogram order
    3. **Poly-phylon genes** — shared across multiple phylons, subclustered
       by Ward's method within each group

    Phylons are ordered by hierarchical clustering dendrogram.

    ## Sorted A_binarized (Phylon-Strain Affinity)

    Strains follow the same three-tier logic:

    1. **No-phylon strains** — not affiliated with any phylon
    2. **Single-phylon strains** — exclusive to one phylon, grouped by phylon
    3. **Multi-phylon strains** — affiliated with multiple phylons, ordered
       by primary affiliation

    Phylons use the same dendrogram order as L_binarized.
    """)
    return


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
def _():
    mo.md("""
    ## Gene Frequency vs Phylon Count

    Linear regression between how frequently a gene appears across genomes
    and the number of phylons it belongs to. Higher-frequency (more conserved)
    genes tend to be active in more phylons, as they are shared across lineages
    rather than being exclusive to one.
    """)
    return


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
    return


@app.cell
def _():
    mo.md("""
    ## Save Outputs
    """)
    return


@app.cell
def _(A_binarized, A_norm, DATA, L_binarized, L_norm, OUT, n_phylons, r2):
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
    return


@app.cell
def _():
    mo.md("""
    ## Per-Phylon Characterization

    Detailed profiles for each phylon — gene counts, exclusive genes,
    functional annotations, and strain demographics — to support
    biological interpretation and phylon naming.
    """)
    return


@app.cell
def _(A_binarized, A_norm, L_binarized, L_norm):
    """Compute per-phylon profile summary table."""
    exclusive_genes_dict = unique_genes_by_phylon(L_binarized)

    _rows = []
    for _p in L_binarized.columns:
        _gene_mask = L_binarized[_p] == 1
        _strain_mask = A_binarized.loc[_p] == 1
        _ng = int(_gene_mask.sum())
        _ns = int(_strain_mask.sum())
        _ne = len(exclusive_genes_dict.get(_p, []))
        _rows.append(
            {
                "phylon": _p,
                "n_genes": _ng,
                "n_strains": _ns,
                "n_exclusive_genes": _ne,
                "pct_exclusive": round(_ne / _ng * 100, 1) if _ng > 0 else 0,
                "mean_L_weight": round(L_norm.loc[_gene_mask, _p].mean(), 4) if _ng > 0 else 0,
                "mean_A_weight": round(A_norm.loc[_p, _strain_mask].mean(), 4) if _ns > 0 else 0,
            }
        )

    phylon_profiles = pd.DataFrame(_rows).sort_values("n_exclusive_genes", ascending=False).reset_index(drop=True)

    mo.output.replace(
        mo.vstack(
            [
                mo.md("### Phylon Profiles (sorted by exclusive gene count)"),
                mo.ui.table(phylon_profiles),
            ]
        )
    )
    return exclusive_genes_dict, phylon_profiles


@app.cell
def _(DATA, L_binarized, SPECIES, exclusive_genes_dict):
    """Map genes to functional annotations via BAKTA/CD-HIT."""
    cluster_products = {}
    phylon_gene_details = None

    try:
        _all_funcs_path = os.path.join(DATA, "processed", "all_functions.csv")
        _bakta_dir = os.path.join(DATA, "processed", "bakta")

        if os.path.exists(_all_funcs_path):
            _funcs = pd.read_csv(_all_funcs_path)
        elif os.path.isdir(_bakta_dir):
            _funcs = collect_functions(DATA, "processed/bakta/")
        else:
            raise FileNotFoundError("No annotation data available")

        _pg_map = get_pg_to_locus_map(DATA, SPECIES)
        _locus_products = _funcs[["locus", "product"]].drop_duplicates()
        _pg_annotated = _pg_map.merge(_locus_products, left_on="gene_id", right_on="locus", how="left")
        cluster_products = (
            _pg_annotated.groupby("cluster")["product"].apply(lambda x: "; ".join(x.dropna().unique())).to_dict()
        )

        # Also build GO mapping
        if "go" in _funcs.columns:
            _locus_go = _funcs[["locus", "go"]].drop_duplicates()
            _pg_go = _pg_map.merge(_locus_go, left_on="gene_id", right_on="locus", how="left")
            _cluster_go = _pg_go.groupby("cluster")["go"].apply(lambda x: "; ".join(x.dropna().unique())).to_dict()
        else:
            _cluster_go = {}

        _rows = []
        for _p in L_binarized.columns:
            _genes = L_binarized.index[L_binarized[_p] == 1].tolist()
            _exclusive = set(exclusive_genes_dict.get(_p, []))
            for _g in _genes:
                _rows.append(
                    {
                        "phylon": _p,
                        "gene": _g,
                        "is_exclusive": _g in _exclusive,
                        "product": cluster_products.get(_g, ""),
                        "go": _cluster_go.get(_g, ""),
                    }
                )
        phylon_gene_details = pd.DataFrame(_rows)

        mo.output.replace(
            mo.md(
                f"Gene annotations loaded: **{len(cluster_products)}** clusters mapped to products\n\n"
                f"Detail table: **{len(phylon_gene_details)}** gene-phylon pairs"
            )
        )
    except Exception as _e:
        mo.output.replace(mo.callout(mo.md(f"Could not load annotations: {_e}"), kind="warn"))

    return cluster_products, phylon_gene_details


@app.cell
def _(A_binarized, metadata):
    """Compute strain demographics (MLST, serovar, host, country) per phylon."""
    _categories = ["mlst", "serovar", "host_name", "isolation_country", "mash_cluster"]
    _available = [c for c in _categories if c in metadata.columns]

    _rows = []
    for _p in A_binarized.index:
        _strain_ids = A_binarized.columns[A_binarized.loc[_p] == 1].tolist()
        _meta_sub = metadata[metadata["genome_id"].isin(_strain_ids)]
        _n = len(_meta_sub)
        if _n == 0:
            continue
        for _cat in _available:
            _counts = _meta_sub[_cat].fillna("Unknown").value_counts()
            for _val, _cnt in _counts.items():
                _rows.append(
                    {
                        "phylon": _p,
                        "category": _cat,
                        "value": str(_val),
                        "count": int(_cnt),
                        "pct": round(_cnt / _n * 100, 1),
                    }
                )

    phylon_strain_demographics = pd.DataFrame(_rows)
    mo.output.replace(
        mo.md(f"Strain demographics: **{len(phylon_strain_demographics)}** entries across {len(_available)} categories")
    )
    return (phylon_strain_demographics,)


@app.cell
def _():
    mo.md("""
    ## Phylon Explorer

    Select a phylon to view its detailed profile: gene counts, top gene
    products, strain MLST/serovar breakdown, and exclusive gene list.
    """)
    return


@app.cell
def _(phylon_profiles):
    """Dropdown to select a phylon for exploration."""
    phylon_dropdown = mo.ui.dropdown(
        options=phylon_profiles["phylon"].tolist(),
        value=phylon_profiles.iloc[0]["phylon"],
        label="Select phylon",
    )
    phylon_dropdown
    return (phylon_dropdown,)


@app.cell
def _(phylon_dropdown, phylon_profiles):
    """Display stat cards for selected phylon."""
    _sel = phylon_dropdown.value
    _row = phylon_profiles[phylon_profiles["phylon"] == _sel].iloc[0]

    mo.vstack(
        [
            mo.md(f"### {_sel}"),
            mo.hstack(
                [
                    mo.stat(value=str(int(_row["n_genes"])), label="Genes"),
                    mo.stat(value=str(int(_row["n_strains"])), label="Strains"),
                    mo.stat(value=str(int(_row["n_exclusive_genes"])), label="Exclusive Genes"),
                    mo.stat(value=f"{_row['pct_exclusive']}%", label="% Exclusive"),
                ],
                justify="center",
                gap="2rem",
            ),
            mo.hstack(
                [
                    mo.stat(value=str(_row["mean_L_weight"]), label="Mean L Weight"),
                    mo.stat(value=str(_row["mean_A_weight"]), label="Mean A Weight"),
                ],
                justify="center",
                gap="2rem",
            ),
        ]
    )
    return


@app.cell
def _(phylon_dropdown, phylon_gene_details, phylon_strain_demographics):
    """Display detailed tabs for selected phylon."""
    _sel = phylon_dropdown.value
    _tabs = {}

    # Top gene products
    if phylon_gene_details is not None:
        _pgd = phylon_gene_details[phylon_gene_details["phylon"] == _sel]
        _products = _pgd["product"].str.split("; ").explode().str.strip()
        _products = _products[(_products != "") & (_products.notna())]
        _prod_counts = _products.value_counts().head(20).reset_index()
        _prod_counts.columns = ["product", "count"]
        _tabs["Top Products"] = mo.ui.table(_prod_counts)

        # Exclusive genes with products
        _excl = _pgd[_pgd["is_exclusive"]].copy()
        if len(_excl) > 0:
            _tabs[f"Exclusive Genes ({len(_excl)})"] = mo.ui.table(_excl[["gene", "product"]].reset_index(drop=True))
        else:
            _tabs["Exclusive Genes (0)"] = mo.md("*No exclusive genes for this phylon.*")

        # All genes
        _tabs[f"All Genes ({len(_pgd)})"] = mo.ui.table(
            _pgd[["gene", "is_exclusive", "product"]].reset_index(drop=True)
        )

    # Strain demographics
    if len(phylon_strain_demographics) > 0:
        _demo = phylon_strain_demographics[phylon_strain_demographics["phylon"] == _sel]
        for _cat in _demo["category"].unique():
            _cat_df = _demo[_demo["category"] == _cat][["value", "count", "pct"]].reset_index(drop=True)
            _tabs[_cat.replace("_", " ").title()] = mo.ui.table(_cat_df)

    if _tabs:
        mo.output.replace(mo.ui.tabs(_tabs))
    else:
        mo.output.replace(mo.md("*No detail data available.*"))
    return


@app.cell
def _():
    mo.md("""
    ## Save Characterization Outputs
    """)
    return


@app.cell
def _(OUT, phylon_gene_details, phylon_profiles, phylon_strain_demographics):
    """Save per-phylon characterization CSVs."""
    _saved = []

    phylon_profiles.to_csv(os.path.join(OUT, "data", "5a_phylon_profiles.csv"), index=False)
    _saved.append(f"- `5a_phylon_profiles.csv`: {len(phylon_profiles)} phylons")

    if phylon_gene_details is not None and len(phylon_gene_details) > 0:
        phylon_gene_details.to_csv(os.path.join(OUT, "data", "5a_phylon_gene_details.csv"), index=False)
        _saved.append(f"- `5a_phylon_gene_details.csv`: {len(phylon_gene_details)} gene-phylon pairs")

    if len(phylon_strain_demographics) > 0:
        phylon_strain_demographics.to_csv(os.path.join(OUT, "data", "5a_phylon_strain_demographics.csv"), index=False)
        _saved.append(f"- `5a_phylon_strain_demographics.csv`: {len(phylon_strain_demographics)} entries")

    mo.output.replace(mo.md("Saved characterization outputs:\n\n" + "\n".join(_saved)))
    return


if __name__ == "__main__":
    app.run()
