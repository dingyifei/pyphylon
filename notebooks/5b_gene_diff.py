import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    import yaml

    from pyphylon.biointerp import collect_functions, get_pg_to_locus_map
    from pyphylon.plotting import generate_phylon_dendrogram
    from pyphylon.plotting_util import find_exclusive_genes

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

    Each dendrogram split is analyzed independently below, showing the two
    branch groups with named phylons, exclusive gene counts, heatmaps, and
    functional annotations.
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
    """Dataset-specific phylon names and split descriptions."""
    PHYLON_NAMES = {
        0: "Phylon 0 \u2014 USA ST-11422 Lineage",
        1: "Phylon 1 \u2014 Finnish ST-677 Clone",
        2: "Phylon 2 \u2014 S. Korean ST-4253 Lineage",
        3: "Phylon 3 \u2014 CC-45/ST-45 Cosmopolitan",
        4: "Phylon 4 \u2014 Mu-like Prophage (MGE)",
        5: "Phylon 5 \u2014 Conjugative T4SS (MGE)",
        6: "Phylon 6 \u2014 Peruvian ST-2993 Lineage",
        7: "Phylon 7 \u2014 CJIE1-like HK97 Phage (MGE)",
        8: "Phylon 8 \u2014 CC-21/ST-50 Poultry Lineage",
        9: "Phylon 9 \u2014 Diverse India/USA Lineage",
        10: "Phylon 10 \u2014 CJIE2 SynExo Prophage (MGE)",
        11: "Phylon 11 \u2014 East Asian Poultry Lineage",
        12: "Phylon 12 \u2014 HS:19/ST-22 GBS-Associated",
    }

    SPLIT_DESCRIPTIONS = {
        "Split 1": (
            "**MGE vs Lineage** \u2014 The root split perfectly separates four "
            "mobile genetic element phylons from nine lineage-defining "
            "chromosomal phylons, reflecting orthogonal evolutionary "
            "dynamics (horizontal vs vertical inheritance)."
        ),
        "Split 2": (
            "**Finnish ST-677 Clone Separates** \u2014 Phylon 1 separates from "
            "all other lineage phylons with 93 ubiquitous exclusive genes "
            "(present in all other lineages but absent from phylon 1), "
            "including flagellin variants, CPS biosynthesis genes, and "
            "autotransporters."
        ),
        "Split 3": (
            "**CC-21 Related vs Remaining Lineages** \u2014 Phylons sharing "
            "CC-21 or CC-21-related ancestry (USA ST-11422, S. Korean "
            "ST-4253, CC-21/ST-50 Poultry) separate from the remaining "
            "lineages."
        ),
        "Split 4": (
            "**India/East Asian vs CC-45/Peru/GBS** \u2014 The Diverse India/USA "
            "and East Asian Poultry lineages separate from the CC-45 "
            "Cosmopolitan, Peruvian ST-2993, and GBS-associated HS:19 "
            "lineages."
        ),
        "Split 5": (
            "**CC-45 Cosmopolitan vs Peru/GBS** \u2014 The extreme generalist "
            "CC-45/ST-45 lineage (pansusceptible, cold-sensitive) separates "
            "from two independent GBS-risk lineages."
        ),
        "Split 6": (
            "**India/USA vs East Asian Poultry** \u2014 The highly diverse "
            "India/USA lineage (40+ STs, T6SS, PEB1 adhesin) separates "
            "from the East Asian Poultry lineage (S. Korea/China focus)."
        ),
        "Split 7": (
            "**USA ST-11422 vs CC-21 Core** \u2014 Geographic/host-specific "
            "divergence within the CC-21 complex: USA-dominant ST-11422 "
            "with transposases vs S. Korean and poultry-associated "
            "sub-lineages."
        ),
        "Split 8": (
            "**Peruvian vs GBS-Associated** \u2014 Two independent GBS-risk "
            "lineages: Phylon 6 (linked to 2019 Peru GBS outbreak) vs "
            "Phylon 12 (HS:19/ST-22 with LOS ganglioside mimicry)."
        ),
        "Split 9": (
            "**S. Korean ST-4253 vs CC-21/ST-50 Poultry** \u2014 S. Korean "
            "lineage with extensive CPS biosynthesis genes separates from "
            "the globally dominant poultry-associated CC-21/ST-50 with "
            "sialyltransferase."
        ),
        "Split 10": (
            "**Mu-like Prophage vs Other MGEs** \u2014 The most widespread MGE "
            "(Mu-like prophage, 28% of strains) separates from the "
            "remaining three mobile genetic elements."
        ),
        "Split 11": (
            "**Conjugative T4SS vs Temperate Phages** \u2014 The conjugative "
            "T4SS element (pVir/pTet-like, 100% exclusive genes, "
            "tetracycline resistance) separates from two temperate phages."
        ),
        "Split 12": (
            "**CJIE1 HK97 vs CJIE2 SynExo** \u2014 Two HK97-family temperate "
            "phages: CJIE1-like with Gam nuclease inhibitor vs CJIE2 with "
            "Bet/YqaJ SynExo recombination module and RusA Holliday "
            "junction resolvase."
        ),
    }

    def name_phylon(col):
        """Look up phylon display name, handling int/str/phylonN column types."""
        if col in PHYLON_NAMES:
            return PHYLON_NAMES[col]
        try:
            return PHYLON_NAMES[int(col)]
        except (ValueError, KeyError):
            pass
        if isinstance(col, str) and col.startswith("phylon"):
            try:
                return PHYLON_NAMES[int(col[6:])]
            except (ValueError, KeyError):
                pass
        return str(col)

    def name_members(membership_str):
        """Convert semicolon-separated membership string to named list."""
        return ", ".join(name_phylon(p) for p in membership_str.split(";"))

    return SPLIT_DESCRIPTIONS, name_members, name_phylon


@app.cell
def _():
    mo.md("""
    ## Load L_binarized

    The binarized L matrix (genes x phylons) assigns each gene to zero or
    more phylons. Values are 0 or 1 from the 3-means thresholding in step 5a.
    """)
    return


@app.cell
def _(DATA):
    """Load L_binarized matrix from 5a outputs."""
    L_binarized = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "L_binarized.csv"), index_col=0)

    mo.output.replace(
        mo.md(
            f"Loaded L_binarized: **{L_binarized.shape[0]}** genes x "
            f"**{L_binarized.shape[1]}** phylons\n\n"
            f"Columns: {list(L_binarized.columns)}"
        )
    )
    return (L_binarized,)


@app.cell
def _():
    mo.md("""
    ## Ward Clustering Heatmap

    Clustermap of L_binarized (genes x phylons) using Ward's minimum-variance
    linkage. Rows and columns are reordered by the dendrogram so that phylons
    with similar gene-membership profiles cluster together.
    """)
    return


@app.cell
def _(FIG, L_binarized, name_phylon):
    """Ward hierarchical clustering heatmap with named phylons."""
    _L_named = L_binarized.rename(columns=lambda c: name_phylon(c))
    g = sns.clustermap(_L_named, method="ward", cmap="hot_r", figsize=(14, 12))

    # Rotate column labels for readability
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8
    )
    g.ax_heatmap.set_yticklabels([])  # too many genes to label rows

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
def _(FIG, L_binarized, name_phylon):
    """Annotated phylon dendrogram with named tick labels."""
    fig_dend, ax_dend = plt.subplots(figsize=(10, 6))
    ax_dend, df_stats, split_gene_sets_labeled = generate_phylon_dendrogram(
        L_binarized,
        text_offset=0,
        labels=["exclusive_genes", "total_split_genes"],
        orientation="left",
        ax=ax_dend,
    )

    # Replace tick labels with descriptive phylon names
    _new_labels = [name_phylon(t.get_text()) for t in ax_dend.get_yticklabels()]
    ax_dend.set_yticklabels(_new_labels, fontsize=7)

    fig_dend.tight_layout()
    fig_dend.savefig(os.path.join(FIG, "5b_phylon_dendrogram.png"), bbox_inches="tight")

    mo.output.replace(mo.vstack([mo.md("**Phylon dendrogram** (annotated with gene counts)"), fig_dend]))
    return df_stats, split_gene_sets_labeled


@app.cell
def _():
    mo.md("""
    ## Interactive Split Explorer
    """)
    return


@app.cell
def _(df_stats):
    """Dropdown to select a dendrogram split for exploration."""
    split_dropdown = mo.ui.dropdown(
        options=df_stats.index.tolist(),
        value=df_stats.index[0],
        label="Select split",
    )
    mo.output.replace(split_dropdown)
    return (split_dropdown,)


@app.cell
def _(df_stats, name_members, split_dropdown, split_gene_sets_labeled):
    _selected = split_dropdown.value
    _row = df_stats.loc[_selected]
    _gene_sets = split_gene_sets_labeled[_selected]

    _stats_row = mo.hstack(
        [
            mo.stat(
                value=str(len(_gene_sets["ubiquitous_exclusive_genes"])),
                label="Ubiquitous Exclusive",
            ),
            mo.stat(value=str(len(_gene_sets["exclusive_genes"])), label="Exclusive"),
            mo.stat(value=str(len(_gene_sets["total_split_genes"])), label="Total Split"),
            mo.stat(
                value=str(len(_gene_sets["total_ubiquitous_genes"])),
                label="Total Ubiquitous",
            ),
        ],
        justify="center",
        gap="2rem",
    )

    _named = name_members(_row["split_membership"])

    mo.vstack(
        [
            mo.md(f"### {_selected}"),
            mo.md(f"**Phylon membership:** {_named}"),
            _stats_row,
        ]
    )
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

    mo.vstack(
        [
            mo.md(f"**Gene lists for {_selected}**"),
            mo.ui.tabs(_tabs),
        ]
    )
    return


@app.cell
def _():
    mo.md("""
    ## Functional Annotations
    """)
    return


@app.cell
def _(DATA, SPECIES):
    """Load BAKTA annotations and build cluster -> product mapping."""
    _all_funcs_path = os.path.join(DATA, "processed", "all_functions.csv")
    _bakta_dir = os.path.join(DATA, "processed", "bakta")

    cluster_products = {}
    _annotation_msg = None

    try:
        if os.path.exists(_all_funcs_path):
            _funcs = pd.read_csv(_all_funcs_path)
        elif os.path.isdir(_bakta_dir):
            _funcs = collect_functions(DATA, "processed/bakta/")
        else:
            _annotation_msg = "BAKTA annotation data not available."
            _funcs = None

        if _funcs is not None:
            _pg_map = get_pg_to_locus_map(DATA, SPECIES)
            _locus_products = _funcs[["locus", "product"]].drop_duplicates()
            _pg_annotated = _pg_map.merge(_locus_products, left_on="gene_id", right_on="locus", how="left")
            cluster_products = (
                _pg_annotated.groupby("cluster")["product"].apply(lambda x: "; ".join(x.dropna().unique())).to_dict()
            )
            _annotation_msg = f"Loaded annotations for **{len(cluster_products)}** gene clusters."
    except Exception as _e:
        _annotation_msg = f"Could not load annotations: {_e}"

    mo.output.replace(mo.md(_annotation_msg if _annotation_msg else "Annotations loaded."))
    return (cluster_products,)


@app.cell
def _(L_binarized, cluster_products, df_stats):
    """Pre-compute exclusive gene analysis for all dendrogram splits."""

    def _find_children(split_name):
        """Find the two direct children of a split in df_stats."""
        parent = set(df_stats.loc[split_name, "split_membership"].split(";"))
        candidates = []
        for label in df_stats.index:
            if label == split_name:
                continue
            members = set(df_stats.loc[label, "split_membership"].split(";"))
            if members < parent:
                candidates.append((label, members))
        # Sort by membership size descending so direct children are found first
        candidates.sort(key=lambda x: len(x[1]), reverse=True)
        for i, (l1, m1) in enumerate(candidates):
            for l2, m2 in candidates[i + 1 :]:
                if m1 | m2 == parent and not (m1 & m2):
                    return l1, sorted(m1), l2, sorted(m2)
        return None, [], None, []

    split_results = {}
    for _split in df_stats.index:
        if not _split.startswith("Split"):
            continue
        _c1_label, _c1_phylons, _c2_label, _c2_phylons = _find_children(_split)
        if _c1_label is None:
            continue

        _g1_only, _g2_only, _shared = find_exclusive_genes(L_binarized, _c1_phylons, _c2_phylons)

        _rows = []
        for _g in _g1_only:
            _rows.append(
                {
                    "split": _split,
                    "branch": _c1_label,
                    "gene": _g,
                    "product": cluster_products.get(_g, ""),
                }
            )
        for _g in _g2_only:
            _rows.append(
                {
                    "split": _split,
                    "branch": _c2_label,
                    "gene": _g,
                    "product": cluster_products.get(_g, ""),
                }
            )

        split_results[_split] = {
            "child1_label": _c1_label,
            "child1_phylons": _c1_phylons,
            "child2_label": _c2_label,
            "child2_phylons": _c2_phylons,
            "group1_only": _g1_only,
            "group2_only": _g2_only,
            "shared": _shared,
            "annotated_rows": _rows,
        }

    mo.output.replace(mo.md(f"Pre-computed analysis for **{len(split_results)}** dendrogram splits."))
    return (split_results,)


@app.cell
def _(SPLIT_DESCRIPTIONS, name_phylon):
    """Reusable function to build a split analysis section."""

    def build_split_section(split_name, split_results, L_binarized, FIG):
        """Build a complete marimo section for one dendrogram split."""
        sr = split_results[split_name]
        c1_label = sr["child1_label"]
        c2_label = sr["child2_label"]
        c1_phylons = sr["child1_phylons"]
        c2_phylons = sr["child2_phylons"]
        g1_only = sr["group1_only"]
        g2_only = sr["group2_only"]
        shared = sr["shared"]

        c1_names = [name_phylon(p) for p in c1_phylons]
        c2_names = [name_phylon(p) for p in c2_phylons]

        split_num = split_name.split()[-1]
        desc = SPLIT_DESCRIPTIONS.get(split_name, "")

        elements = []

        # Header
        elements.append(mo.md(f"### {split_name}"))

        # Description
        if desc:
            elements.append(mo.md(desc))

        # Branch callouts with named phylons
        elements.append(
            mo.hstack(
                [
                    mo.callout(
                        mo.md(f"**{c1_label}**\n\n" + "\n".join(f"- {n}" for n in c1_names)),
                        kind="info",
                    ),
                    mo.callout(
                        mo.md(f"**{c2_label}**\n\n" + "\n".join(f"- {n}" for n in c2_names)),
                        kind="info",
                    ),
                ],
                justify="center",
                gap="2rem",
            )
        )

        # Exclusive gene stats
        elements.append(
            mo.hstack(
                [
                    mo.stat(value=str(len(g1_only)), label=f"{c1_label} exclusive"),
                    mo.stat(value=str(len(g2_only)), label=f"{c2_label} exclusive"),
                    mo.stat(value=str(len(shared)), label="Shared"),
                ],
                justify="center",
                gap="2rem",
            )
        )

        # Heatmap of exclusive genes
        _gene_rows = list(g1_only) + list(g2_only)
        _cols = c1_phylons + c2_phylons
        if len(_gene_rows) > 0:
            _heat_df = L_binarized.loc[_gene_rows, _cols].copy()
            _heat_df.columns = [name_phylon(c) for c in _heat_df.columns]
            _fig, _ax = plt.subplots(
                figsize=(
                    max(6, len(_cols) * 1.2),
                    max(4, len(_gene_rows) * 0.08),
                )
            )
            sns.heatmap(
                _heat_df,
                cmap="hot_r",
                cbar_kws={"label": "Present"},
                yticklabels=len(_gene_rows) <= 100,
                ax=_ax,
            )
            _ax.axvline(x=len(c1_phylons), color="red", linewidth=2)
            _ax.axhline(y=len(g1_only), color="blue", linewidth=1, linestyle="--")
            _ax.set_xticklabels(_ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
            _ax.set_title(f"Exclusive genes: {c1_label} vs {c2_label}")
            _fig.tight_layout()
            _fig.savefig(
                os.path.join(FIG, f"5b_split_{split_num}_heatmap.png"),
                bbox_inches="tight",
            )
            elements.append(_fig)

        # Top-15 product descriptions per branch
        annotated = sr["annotated_rows"]
        if annotated:
            _df_ann = pd.DataFrame(annotated)
            _tabs = {}
            for _branch in [c1_label, c2_label]:
                _branch_df = _df_ann[_df_ann["branch"] == _branch]
                _products = _branch_df["product"].str.split("; ").explode().str.strip()
                _products = _products[_products != ""].value_counts().head(15)
                if len(_products) > 0:
                    _pdf = _products.reset_index()
                    _pdf.columns = ["product", "count"]
                    _tabs[f"{_branch} ({len(_branch_df)} genes)"] = mo.ui.table(_pdf)
                else:
                    _tabs[f"{_branch} (0)"] = mo.md("*No annotated products.*")
            elements.append(mo.md("**Top-15 product descriptions per branch**"))
            elements.append(mo.ui.tabs(_tabs))

        return mo.vstack(elements)

    return (build_split_section,)


@app.cell
def _():
    mo.md("""
    ## Dendrogram Split Analysis

    Each dendrogram split is presented below as an independent section.
    For every split the two child branches are shown with their named
    phylons, exclusive gene counts, a heatmap of differentiating genes,
    and the top functional annotations per branch.
    """)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 1: MGE cluster vs Lineage cluster (root split)."""
    build_split_section("Split 1", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 2: Finnish ST-677 Clone separates."""
    build_split_section("Split 2", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 3: CC-21 related vs remaining lineages."""
    build_split_section("Split 3", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 4: India/East Asian vs CC-45/Peru/GBS."""
    build_split_section("Split 4", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 5: CC-45 Cosmopolitan vs Peru/GBS."""
    build_split_section("Split 5", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 6: India/USA vs East Asian Poultry."""
    build_split_section("Split 6", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 7: USA ST-11422 vs CC-21 core."""
    build_split_section("Split 7", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 8: Peruvian vs GBS-associated."""
    build_split_section("Split 8", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 9: S. Korean ST-4253 vs CC-21/ST-50 Poultry."""
    build_split_section("Split 9", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 10: Mu-like Prophage vs other MGEs."""
    build_split_section("Split 10", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 11: Conjugative T4SS vs temperate phages."""
    build_split_section("Split 11", split_results, L_binarized, FIG)
    return


@app.cell
def _(FIG, L_binarized, build_split_section, split_results):
    """Split 12: CJIE1 HK97 vs CJIE2 SynExo."""
    build_split_section("Split 12", split_results, L_binarized, FIG)
    return


@app.cell
def _():
    mo.md("""
    ## Save Results
    """)
    return


@app.cell
def _(L_binarized, OUT, df_stats, split_results):
    """Save dendrogram statistics and comprehensive split gene annotations."""
    df_stats.to_csv(os.path.join(OUT, "data", "5b_gene_diff_stats.csv"))

    _all_rows = []
    for _sr in split_results.values():
        _all_rows.extend(_sr["annotated_rows"])

    _extra_lines = []
    if _all_rows:
        _df_all = pd.DataFrame(_all_rows)
        _df_all.to_csv(os.path.join(OUT, "data", "5b_all_split_genes.csv"), index=False)
        _extra_lines.append(
            f"- `5b_all_split_genes.csv`: {len(_df_all)} annotated genes across {len(split_results)} splits"
        )

    _heatmap_count = sum(1 for _sr in split_results.values() if len(_sr["group1_only"]) + len(_sr["group2_only"]) > 0)

    mo.output.replace(
        mo.md(
            f"Saved outputs:\n\n"
            f"- `5b_gene_diff_stats.csv`: {df_stats.shape[0]} rows "
            f"(splits + leaves)\n"
            f"- `5b_clustermap.png`: Ward clustering heatmap\n"
            f"- `5b_phylon_dendrogram.png`: annotated dendrogram\n"
            f"- `5b_split_*_heatmap.png`: {_heatmap_count} per-split heatmaps\n"
            + ("\n".join(_extra_lines) + "\n" if _extra_lines else "")
            + f"\nL_binarized: {L_binarized.shape[0]} genes x "
            f"{L_binarized.shape[1]} phylons"
        )
    )
    return


if __name__ == "__main__":
    app.run()
