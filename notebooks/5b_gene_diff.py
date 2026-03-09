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
    from pyphylon.plotting_util import find_exclusive_genes
    from pyphylon.biointerp import collect_functions, get_pg_to_locus_map

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    mo.md("""
    # 5b: Gene Differentiation

    Hierarchical clustering (**Ward's method**) of the binarized L matrix
    reveals relationships between phylons based on shared gene content.
    The companion dendrogram annotates each split with **exclusive gene
    counts** and **total split gene counts**, quantifying how gene
    repertoires diverge across phylon lineages.
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

    DATA = CONFIG["DATA_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]
    SPECIES = CONFIG["PG_NAME"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    return DATA, FIG, OUT, SPECIES


@app.cell
def _():
    mo.md("""
    ## Load L_binarized

    The binarized L matrix (genes × phylons) assigns each gene to zero or
    more phylons. Values are 0 or 1 from the 3-means thresholding in step 5a.
    """)
    return


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
def _():
    mo.md("""
    ## Ward Clustering Heatmap

    Clustermap of L_binarized (genes × phylons) using Ward's minimum-variance
    linkage. Rows and columns are reordered by the dendrogram so that phylons
    with similar gene-membership profiles cluster together.
    """)
    return


@app.cell
def _(FIG, L_binarized):
    """Ward hierarchical clustering heatmap."""
    g = sns.clustermap(L_binarized, method="ward", cmap="hot_r")
    g.savefig(os.path.join(FIG, "5b_clustermap.png"), bbox_inches="tight")
    mo.output.replace(g.fig)
    return


@app.cell
def _():
    mo.md("""
    ## Phylon Dendrogram

    Annotated dendrogram showing gene differentiation at each split node.
    Labels indicate **exclusive genes** (genes unique to one branch) and
    **total split genes** (all genes distinguishing the two branches).
    """)
    return


@app.cell
def _(FIG, L_binarized):
    """Annotated phylon dendrogram with gene counts per split."""
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
def _(df_stats):
    """Dropdown to select a dendrogram split for exploration."""
    split_dropdown = mo.ui.dropdown(
        options=df_stats.index.tolist(),
        value=df_stats.index[0],
        label="Select split",
    )
    split_dropdown
    return (split_dropdown,)


@app.cell
def _(df_stats, split_dropdown, split_gene_sets_labeled):
    _selected = split_dropdown.value
    _row = df_stats.loc[_selected]
    _gene_sets = split_gene_sets_labeled[_selected]

    _stats_row = mo.hstack(
        [
            mo.stat(value=str(len(_gene_sets["ubiquitous_exclusive_genes"])), label="Ubiquitous Exclusive"),
            mo.stat(value=str(len(_gene_sets["exclusive_genes"])), label="Exclusive"),
            mo.stat(value=str(len(_gene_sets["total_split_genes"])), label="Total Split"),
            mo.stat(value=str(len(_gene_sets["total_ubiquitous_genes"])), label="Total Ubiquitous"),
        ],
        justify="center",
        gap="2rem",
    )

    _membership = _row["split_membership"]

    mo.vstack([
        mo.md(f"### {_selected}"),
        mo.md(f"**Phylon membership:** {_membership}"),
        _stats_row,
    ])
    return


@app.cell
def _(split_dropdown, split_gene_sets_labeled):
    _selected = split_dropdown.value
    _gene_sets = split_gene_sets_labeled[_selected]

    _category_labels = {
        "ubiquitous_exclusive_genes": "Ubiquitous Exclusive",
        "exclusive_genes": "Exclusive",
        "total_split_genes": "Total Split",
        "total_ubiquitous_genes": "Total Ubiquitous",
    }

    _tabs = {}
    for _key, _label in _category_labels.items():
        _genes = _gene_sets[_key]
        _df = pd.DataFrame({"gene": list(_genes)})
        if len(_df) > 0:
            _tabs[f"{_label} ({len(_genes)})"] = mo.ui.table(_df)
        else:
            _tabs[f"{_label} (0)"] = mo.md("*No genes in this category.*")

    mo.vstack([
        mo.md(f"**Gene lists for {_selected}**"),
        mo.ui.tabs(_tabs),
    ])
    return


@app.cell
def _():
    mo.md("""
    ## First Split: Branch-vs-Branch Comparison

    The interactive explorer above shows gene sets **per cluster** (one side
    vs all other phylons globally). This section compares the two branches
    of the **first dendrogram split** directly, identifying genes that
    distinguish one branch from the other.
    """)
    return


@app.cell
def _(L_binarized, df_stats):
    """Extract first-split children and run find_exclusive_genes."""
    # Split 1 is the root — its membership covers all phylons.
    _root_phylons = set(df_stats.loc["Split 1", "split_membership"].split(";"))

    # Walk subsequent df_stats rows to find the two direct children
    # (the first two entries whose memberships are non-overlapping subsets of the root).
    _child1_phylons = None
    _child2_phylons = None
    _child1_label = None
    _child2_label = None

    for _label in df_stats.index[1:]:
        _members = set(df_stats.loc[_label, "split_membership"].split(";"))
        if not _members.issubset(_root_phylons):
            continue
        if _child1_phylons is None:
            _child1_phylons = _members
            _child1_label = _label
        elif _members.issubset(_root_phylons - _child1_phylons):
            _child2_phylons = _members
            _child2_label = _label
            break

    child1_phylons = sorted(_child1_phylons)
    child2_phylons = sorted(_child2_phylons)
    child1_label = _child1_label
    child2_label = _child2_label

    group1_only, group2_only, shared_genes = find_exclusive_genes(
        L_binarized, child1_phylons, child2_phylons
    )

    mo.output.replace(
        mo.vstack([
            mo.md("### First Split Children"),
            mo.hstack([
                mo.callout(
                    mo.md(f"**{child1_label}**\n\n{', '.join(child1_phylons)}"),
                    kind="info",
                ),
                mo.callout(
                    mo.md(f"**{child2_label}**\n\n{', '.join(child2_phylons)}"),
                    kind="info",
                ),
            ], justify="center", gap="2rem"),
            mo.hstack([
                mo.stat(value=str(len(group1_only)), label=f"{child1_label} exclusive"),
                mo.stat(value=str(len(group2_only)), label=f"{child2_label} exclusive"),
                mo.stat(value=str(len(shared_genes)), label="Shared"),
            ], justify="center", gap="2rem"),
        ])
    )
    return (
        child1_label,
        child1_phylons,
        child2_label,
        child2_phylons,
        group1_only,
        group2_only,
    )


@app.cell
def _(
    FIG,
    L_binarized,
    child1_label,
    child1_phylons,
    child2_label,
    child2_phylons,
    group1_only,
    group2_only,
):
    """Heatmap of exclusive genes for the first split."""
    _rows = list(group1_only) + list(group2_only)
    _cols = child1_phylons + child2_phylons

    if len(_rows) == 0:
        mo.output.replace(mo.md("*No exclusive genes found for the first split.*"))
    else:
        _heat_df = L_binarized.loc[_rows, _cols]
        _fig, _ax = plt.subplots(figsize=(max(6, len(_cols) * 0.8), max(4, len(_rows) * 0.08)))
        sns.heatmap(_heat_df, cmap="hot_r", cbar_kws={"label": "Present"}, ax=_ax)
        _ax.axvline(x=len(child1_phylons), color="red", linewidth=2)
        _ax.axhline(y=len(group1_only), color="blue", linewidth=1, linestyle="--")
        _ax.set_title(f"Exclusive genes: {child1_label} vs {child2_label}")
        _fig.tight_layout()
        _fig.savefig(os.path.join(FIG, "5b_first_split_heatmap.png"), bbox_inches="tight")
        mo.output.replace(_fig)
    return


@app.cell
def _(DATA, SPECIES, child1_label, child2_label, group1_only, group2_only):
    """Annotate exclusive genes with BAKTA functional descriptions."""
    _all_funcs_path = os.path.join(DATA, "processed", "all_functions.csv")
    _bakta_dir = os.path.join(DATA, "processed", "bakta")

    df_annotated = None
    _skip_msg = None

    try:
        if os.path.exists(_all_funcs_path):
            _funcs = pd.read_csv(_all_funcs_path)
        elif os.path.isdir(_bakta_dir):
            _funcs = collect_functions(DATA, "processed/bakta/")
        else:
            _skip_msg = "BAKTA annotation data not available — skipping functional annotation."
            _funcs = None

        if _funcs is not None:
            _pg_map = get_pg_to_locus_map(DATA, SPECIES)
            # Map cluster → product via locus join
            _locus_products = _funcs[["locus", "product"]].drop_duplicates()
            _pg_annotated = _pg_map.merge(_locus_products, left_on="gene_id", right_on="locus", how="left")
            _cluster_products = (
                _pg_annotated.groupby("cluster")["product"]
                .apply(lambda x: "; ".join(x.dropna().unique()))
                .to_dict()
            )

            _rows = []
            for _gene in group1_only:
                _rows.append({"gene": _gene, "branch": child1_label, "product": _cluster_products.get(_gene, "")})
            for _gene in group2_only:
                _rows.append({"gene": _gene, "branch": child2_label, "product": _cluster_products.get(_gene, "")})
            df_annotated = pd.DataFrame(_rows)
    except Exception as _e:
        _skip_msg = f"Could not load annotations: {_e}"

    if _skip_msg:
        mo.output.replace(mo.callout(mo.md(_skip_msg), kind="warn"))
    elif df_annotated is not None and len(df_annotated) > 0:
        mo.output.replace(
            mo.vstack([
                mo.md(f"**Annotated exclusive genes** ({len(df_annotated)} genes)"),
                mo.ui.table(df_annotated),
            ])
        )
    return (df_annotated,)


@app.cell
def _(child1_label, child2_label, df_annotated):
    """Top-15 product descriptions per branch."""
    if df_annotated is None or len(df_annotated) == 0 or "product" not in df_annotated.columns:
        mo.output.replace(mo.md("*No annotation data to summarize.*"))
    else:
        _tabs = {}
        for _branch in [child1_label, child2_label]:
            _branch_df = df_annotated[df_annotated["branch"] == _branch]
            # Explode semicolon-joined products and count
            _products = _branch_df["product"].str.split("; ").explode().str.strip()
            _products = _products[_products != ""].value_counts().head(15)
            if len(_products) > 0:
                _pdf = _products.reset_index()
                _pdf.columns = ["product", "count"]
                _tabs[f"{_branch} ({len(_branch_df)} genes)"] = mo.ui.table(_pdf)
            else:
                _tabs[f"{_branch} (0)"] = mo.md("*No annotated products.*")
        mo.output.replace(
            mo.vstack([
                mo.md("**Top-15 product descriptions per branch**"),
                mo.ui.tabs(_tabs),
            ])
        )
    return


@app.cell
def _():
    mo.md("""
    ## Save Results
    """)
    return


@app.cell
def _(L_binarized, OUT, df_annotated, df_stats):
    """Save dendrogram statistics and first-split gene annotations."""
    df_stats.to_csv(os.path.join(OUT, "data", "5b_gene_diff_stats.csv"))

    _extra_lines = []
    if df_annotated is not None and len(df_annotated) > 0:
        df_annotated.to_csv(os.path.join(OUT, "data", "5b_first_split_genes.csv"), index=False)
        _extra_lines.append(f"- `5b_first_split_genes.csv`: {len(df_annotated)} annotated exclusive genes")

    mo.output.replace(
        mo.md(
            f"Saved outputs:\n\n"
            f"- `5b_gene_diff_stats.csv`: {df_stats.shape[0]} rows (splits + leaves)\n"
            f"- `5b_clustermap.png`: Ward clustering heatmap\n"
            f"- `5b_phylon_dendrogram.png`: annotated dendrogram\n"
            f"- `5b_first_split_heatmap.png`: first-split exclusive genes heatmap\n"
            + ("\n".join(_extra_lines) + "\n" if _extra_lines else "")
            + f"\nL_binarized: {L_binarized.shape[0]} genes x {L_binarized.shape[1]} phylons"
        )
    )
    return


if __name__ == "__main__":
    app.run()
