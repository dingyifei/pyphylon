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
    from scipy import stats

    from pyphylon.biointerp import (
        calc_all_phylon_go_enrichments,
        collect_functions,
        explode_go_annos,
        gen_phylon_wordcloud,
        get_go_mapping,
        get_pg_to_locus_map,
        load_cogclassifier_results,
        load_eggnog_annotations,
    )

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    mo.md("""
    # 5c: Functional Enrichments

    GO term enrichment analysis identifies biological functions that are
    significantly over-represented in each phylon's gene set. For every
    phylon-GO term pair a **hypergeometric test** compares the observed
    overlap to the expected overlap under random sampling; results are
    filtered at **p < 0.05**.

    The generic category *SO:0001217* ("protein-coding gene") is excluded
    because it carries no discriminating information.
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
    ## Load Inputs
    """)
    return


@app.cell
def _(DATA):
    """Load L_binarized matrix and rename columns to 1-indexed phylon names."""
    L_BIN = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "L_binarized.csv"), index_col=0)
    L_BIN.columns = [f"phylon{i}" for i in range(1, L_BIN.shape[1] + 1)]

    mo.output.replace(mo.md(f"Loaded L_binarized: **{L_BIN.shape[0]}** genes x **{L_BIN.shape[1]}** phylons"))
    return (L_BIN,)


@app.cell
def _(DATA):
    """Collect gene functions from BAKTA annotations (existence-check pattern)."""
    all_functions_path = os.path.join(DATA, "processed", "all_functions.csv")

    if not os.path.exists(all_functions_path):
        all_functions_raw = collect_functions(DATA, "processed/bakta/")
        all_functions_raw.to_csv(all_functions_path)

    all_functions = pd.read_csv(all_functions_path, index_col=0)

    mo.output.replace(mo.md(f"Loaded **{len(all_functions)}** function annotations from BAKTA"))
    return (all_functions,)


@app.cell
def _():
    mo.md("""
    ## Build GO Mapping

    Each CD-HIT gene cluster is linked to its BAKTA-annotated product
    description and GO terms via the pangenome-to-locus map. GO terms
    assigned to fewer than 4 clusters are dropped to reduce noise.
    """)
    return


@app.cell
def _(DATA, SPECIES, all_functions):
    """Build cluster → function mapping and filter GO terms."""
    pg2locus_map = get_pg_to_locus_map(DATA, SPECIES)
    functions2genes = pd.merge(all_functions, pg2locus_map, left_on="locus", right_on="gene_id")

    cluster_functions = functions2genes.groupby("cluster").first().reset_index()[["cluster", "product", "go"]]
    cluster_to_go_functions = explode_go_annos(cluster_functions)

    go_functions_count = cluster_to_go_functions.groupby("go").count()
    go_functions = go_functions_count[go_functions_count["cluster"] > 3].sort_values("cluster", ascending=False)

    mo.output.replace(
        mo.md(
            f"Cluster → function mapping:\n\n"
            f"- **{len(pg2locus_map)}** cluster-to-locus mappings\n"
            f"- **{len(functions2genes)}** merged function-gene records\n"
            f"- **{len(cluster_to_go_functions)}** cluster-GO pairs (exploded)\n"
            f"- **{len(go_functions)}** GO terms with >3 clusters (used for enrichment)"
        )
    )
    return cluster_to_go_functions, functions2genes, go_functions


@app.cell
def _():
    mo.md("""
    ## Compute Enrichments
    """)
    return


@app.cell
def _(L_BIN, OUT, cluster_to_go_functions, functions2genes, go_functions):
    """Compute GO enrichments per phylon (cached), filter, merge GO names, save."""
    with mo.persistent_cache("5c_enrichments"):
        phylon_go_enrichments_raw = calc_all_phylon_go_enrichments(
            L_BIN, functions2genes, cluster_to_go_functions, go_functions, phylon_contribution_cutoff=0.5
        )

    # Filter to significant enrichments
    phylon_go_enrichments = phylon_go_enrichments_raw[phylon_go_enrichments_raw["p_value"] < 0.05].copy()

    # Merge human-readable GO names
    go_mapping = get_go_mapping()
    phylon_go_enrichments = pd.merge(
        phylon_go_enrichments, go_mapping, left_on="function", right_index=True, how="left"
    )
    missing_go = phylon_go_enrichments[phylon_go_enrichments["name"].isnull()].index
    phylon_go_enrichments.loc[missing_go, "name"] = phylon_go_enrichments.loc[missing_go, "function"]

    # Filter out SO:0001217 (generic "protein encoding gene" category)
    phylon_go_enrichments = phylon_go_enrichments[phylon_go_enrichments["function"] != "SO:0001217"]

    # Save enrichment results
    phylon_go_enrichments.to_csv(os.path.join(OUT, "data", "5c_phylon_go_enrichment.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"GO enrichments (p < 0.05):\n\n"
            f"- **{len(phylon_go_enrichments_raw)}** total enrichment tests\n"
            f"- **{len(phylon_go_enrichments)}** significant (p < 0.05, excl. SO:0001217)\n"
            f"- **{phylon_go_enrichments['phylon'].nunique()}** phylons with enrichments\n"
            f"- **{phylon_go_enrichments['function'].nunique()}** unique GO terms enriched"
        )
    )
    return (phylon_go_enrichments,)


@app.cell
def _():
    mo.md("""
    ## Top Enrichments
    """)
    return


@app.cell
def _(phylon_go_enrichments):
    """Interactive table of significant GO enrichments."""
    display_cols = ["phylon", "function", "name", "p_value", "overlap", "genes in phylon", "genes w function"]
    available_cols = [c for c in display_cols if c in phylon_go_enrichments.columns]
    mo.ui.table(
        phylon_go_enrichments[available_cols].sort_values("p_value"),
        label="Significant GO enrichments (p < 0.05)",
    )
    return


@app.cell
def _():
    mo.md("""
    ## Enrichment Heatmap

    Clustered heatmap of enrichment p-values (phylons x GO terms).
    Darker colours indicate stronger enrichment (lower p-values);
    missing values are filled with 0.05 (the significance threshold).
    """)
    return


@app.cell
def _(FIG, phylon_go_enrichments):
    """Generate enrichment heatmap (phylon × GO term, colored by p-value)."""
    phylon_go_enrichments_mat = pd.pivot_table(
        phylon_go_enrichments, index="phylon", columns="function", values="p_value"
    )
    g = sns.clustermap(phylon_go_enrichments_mat.fillna(0.05), cmap="rocket_r")
    g.savefig(os.path.join(FIG, "5c_enrichment_heatmap.png"), bbox_inches="tight")
    plt.close(g.fig)

    mo.output.replace(mo.md("Saved `5c_enrichment_heatmap.png`"))
    return


@app.cell
def _():
    mo.md("""
    ## Phylon Wordclouds

    Word clouds of gene product descriptions for each phylon, weighted
    by phylon contribution. Common uninformative terms (hypothetical,
    uncharacterised, etc.) are filtered out by `gen_phylon_wordcloud`.
    """)
    return


@app.cell
def _(FIG, L_BIN, functions2genes, phylon_go_enrichments):
    """Generate wordclouds per phylon and display in tabs."""
    wordcloud_phylons = sorted(phylon_go_enrichments["phylon"].unique())
    wc_tabs = {}
    for _phylon in wordcloud_phylons:
        wc_path = os.path.join(FIG, f"5c_wordcloud_{_phylon}.png")
        gen_phylon_wordcloud(L_BIN, functions2genes, _phylon, cutoff=0, save=True, filename=wc_path)
        wc_tabs[_phylon] = mo.image(src=wc_path)

    mo.ui.tabs(wc_tabs)
    return


@app.cell
def _():
    mo.md("""
    ## COG & KEGG Enrichment (eggNOG-mapper + COGclassifier)

    COG functional categories from COGclassifier (official NCBI COG database
    via RPS-BLAST) and KEGG pathway annotations from eggNOG-mapper provide
    complementary functional views of each phylon's gene content.
    """)
    return


@app.cell
def _(DATA, L_BIN, SPECIES):
    """Load eggNOG-mapper and COGclassifier annotations."""
    eggnog_annos = None
    cog_annos = None

    try:
        eggnog_annos = load_eggnog_annotations(DATA, SPECIES)
        _eggnog_in_L = eggnog_annos.index.isin(L_BIN.index).sum()
        _eggnog_cog = (eggnog_annos["cog_category"] != "").sum()
        _eggnog_kegg = (eggnog_annos["kegg_ko"] != "").sum()
    except FileNotFoundError:
        _eggnog_in_L = 0
        _eggnog_cog = 0
        _eggnog_kegg = 0

    try:
        cog_annos = load_cogclassifier_results(DATA, SPECIES)
        _cog_in_L = cog_annos.index.isin(L_BIN.index).sum()
        _cog_total = len(cog_annos)
    except FileNotFoundError:
        _cog_in_L = 0
        _cog_total = 0

    _msgs = []
    if eggnog_annos is not None:
        _msgs.append(
            f"- **eggNOG-mapper:** {len(eggnog_annos)} clusters annotated "
            f"({_eggnog_in_L} in L_binarized), "
            f"{_eggnog_cog} with COG, {_eggnog_kegg} with KEGG KO"
        )
    else:
        _msgs.append("- **eggNOG-mapper:** not available (run `snakemake eggnog_annotate --sdm conda`)")

    if cog_annos is not None:
        _msgs.append(
            f"- **COGclassifier:** {_cog_total} clusters annotated "
            f"({_cog_in_L} in L_binarized)"
        )
    else:
        _msgs.append("- **COGclassifier:** not available (run `snakemake cogclassifier_run --sdm conda`)")

    mo.output.replace(mo.md("Annotation sources:\n\n" + "\n".join(_msgs)))
    return cog_annos, eggnog_annos


@app.cell
def _(FIG, L_BIN, OUT, cog_annos, eggnog_annos):
    """COG category distribution and enrichment per phylon."""
    # COG category full names
    _cog_names = {
        "A": "RNA processing", "B": "Chromatin", "C": "Energy",
        "D": "Cell cycle", "E": "Amino acid", "F": "Nucleotide",
        "G": "Carbohydrate", "H": "Coenzyme", "I": "Lipid",
        "J": "Translation", "K": "Transcription", "L": "Replication",
        "M": "Cell wall", "N": "Cell motility", "O": "Post-translational",
        "P": "Inorganic ion", "Q": "Secondary metabolites",
        "S": "Function unknown", "T": "Signal transduction",
        "U": "Secretion", "V": "Defense", "W": "Extracellular",
        "X": "Mobilome", "Y": "Nuclear", "Z": "Cytoskeleton",
    }

    # Prefer COGclassifier (official NCBI COGs); fall back to eggNOG
    if cog_annos is not None:
        _cog_source = cog_annos[["cog_category"]].copy()
        _cog_label = "COGclassifier"
    elif eggnog_annos is not None:
        _cog_source = eggnog_annos[["cog_category"]].copy()
        _cog_label = "eggNOG-mapper"
    else:
        mo.output.replace(mo.callout(mo.md("No COG annotations available."), kind="warn"))
        phylon_cog_enrichment = pd.DataFrame()
        return (phylon_cog_enrichment,)

    # Restrict to genes in L_binarized
    _cog_source = _cog_source[_cog_source.index.isin(L_BIN.index)]
    _cog_source = _cog_source[_cog_source["cog_category"] != ""]

    # Explode multi-letter COG categories (e.g. "KL" -> ["K", "L"])
    _cog_exploded = _cog_source.copy()
    _cog_exploded["cog_category"] = _cog_exploded["cog_category"].apply(list)
    _cog_exploded = _cog_exploded.explode("cog_category")

    # Per-phylon COG counts
    _rows = []
    for _p in L_BIN.columns:
        _genes = L_BIN.index[L_BIN[_p] == 1]
        _cog_in_phylon = _cog_exploded[_cog_exploded.index.isin(_genes)]
        _vc = _cog_in_phylon["cog_category"].value_counts()
        for _cat, _cnt in _vc.items():
            _rows.append({"phylon": _p, "cog_category": _cat, "count": _cnt})

    _cog_phylon_df = pd.DataFrame(_rows)
    if _cog_phylon_df.empty:
        mo.output.replace(mo.callout(mo.md("No COG categories mapped to phylon genes."), kind="warn"))
        phylon_cog_enrichment = pd.DataFrame()
        return (phylon_cog_enrichment,)

    # Pivot for heatmap
    _cog_matrix = _cog_phylon_df.pivot_table(
        index="phylon", columns="cog_category", values="count", fill_value=0
    )
    # Sort columns by total count
    _cog_matrix = _cog_matrix[_cog_matrix.sum().sort_values(ascending=False).index]

    # Hypergeometric enrichment test per phylon-COG pair
    _N = len(_cog_exploded)  # total cluster-COG assignments
    _enrich_rows = []
    for _p in L_BIN.columns:
        _genes = set(L_BIN.index[L_BIN[_p] == 1])
        _n = len(_cog_exploded[_cog_exploded.index.isin(_genes)])  # genes in phylon with any COG
        for _cat in _cog_matrix.columns:
            _K = int((_cog_exploded["cog_category"] == _cat).sum())  # total with this COG
            _k = int(_cog_matrix.loc[_p, _cat]) if _p in _cog_matrix.index else 0
            _pval = stats.hypergeom.sf(_k - 1, _N, _K, _n) if _k > 0 else 1.0
            _enrich_rows.append({
                "phylon": _p,
                "cog_category": _cat,
                "cog_name": _cog_names.get(_cat, _cat),
                "count": _k,
                "total_in_category": _K,
                "genes_in_phylon": _n,
                "p_value": _pval,
            })
    phylon_cog_enrichment = pd.DataFrame(_enrich_rows)
    phylon_cog_enrichment.to_csv(os.path.join(OUT, "data", "5c_phylon_cog_enrichment.csv"), index=False)

    # --- Plots ---
    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [1, 1.2]})

    # Heatmap
    sns.heatmap(
        _cog_matrix,
        annot=True,
        fmt=".0f",
        cmap="YlOrRd",
        ax=_ax1,
        cbar_kws={"label": "Gene count"},
        linewidths=0.5,
    )
    _ax1.set_title(f"COG Categories by Phylon ({_cog_label})")
    _ax1.set_ylabel("")
    _ax1.set_xlabel("COG category")

    # Stacked bar
    _pct_matrix = _cog_matrix.div(_cog_matrix.sum(axis=1), axis=0) * 100
    _palette = sns.color_palette("tab20", n_colors=len(_pct_matrix.columns))
    _left = np.zeros(len(_pct_matrix))
    for _i, _col in enumerate(_pct_matrix.columns):
        _ax2.barh(
            _pct_matrix.index, _pct_matrix[_col], left=_left,
            color=_palette[_i], label=f"{_col} ({_cog_names.get(_col, '')})",
            edgecolor="white", linewidth=0.3,
        )
        _left += _pct_matrix[_col].values
    _ax2.set_xlabel("% of phylon COG assignments")
    _ax2.set_title("COG Composition by Phylon")
    _ax2.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, ncol=1)
    _ax2.invert_yaxis()

    plt.tight_layout()
    _fig.savefig(os.path.join(FIG, "5c_cog_heatmap.png"), bbox_inches="tight")
    mo.output.replace(_fig)
    return (phylon_cog_enrichment,)


@app.cell
def _(FIG, L_BIN, OUT, eggnog_annos):
    """KEGG pathway enrichment per phylon."""
    if eggnog_annos is None or "kegg_pathway" not in eggnog_annos.columns:
        mo.output.replace(mo.callout(mo.md("No KEGG pathway annotations available."), kind="warn"))
        phylon_kegg_enrichment = pd.DataFrame()
        return (phylon_kegg_enrichment,)

    # Restrict to L_binarized genes with KEGG pathways
    _kegg = eggnog_annos[eggnog_annos.index.isin(L_BIN.index)][["kegg_pathway"]].copy()
    _kegg = _kegg[_kegg["kegg_pathway"] != ""]

    # Explode comma-separated pathways
    _kegg["kegg_pathway"] = _kegg["kegg_pathway"].str.split(",")
    _kegg_exploded = _kegg.explode("kegg_pathway")
    _kegg_exploded["kegg_pathway"] = _kegg_exploded["kegg_pathway"].str.strip()
    _kegg_exploded = _kegg_exploded[_kegg_exploded["kegg_pathway"] != ""]

    # Filter to pathways with >3 cluster hits
    _pw_counts = _kegg_exploded["kegg_pathway"].value_counts()
    _valid_pws = _pw_counts[_pw_counts > 3].index
    _kegg_exploded = _kegg_exploded[_kegg_exploded["kegg_pathway"].isin(_valid_pws)]

    if _kegg_exploded.empty:
        mo.output.replace(mo.callout(mo.md("No KEGG pathways with >3 hits."), kind="warn"))
        phylon_kegg_enrichment = pd.DataFrame()
        return (phylon_kegg_enrichment,)

    # Hypergeometric enrichment
    _N = len(_kegg_exploded)
    _enrich_rows = []
    for _p in L_BIN.columns:
        _genes = set(L_BIN.index[L_BIN[_p] == 1])
        _n = len(_kegg_exploded[_kegg_exploded.index.isin(_genes)])
        for _pw in _valid_pws:
            _K = int((_kegg_exploded["kegg_pathway"] == _pw).sum())
            _k = int(_kegg_exploded[_kegg_exploded.index.isin(_genes) & (_kegg_exploded["kegg_pathway"] == _pw)].shape[0])
            _pval = stats.hypergeom.sf(_k - 1, _N, _K, _n) if _k > 0 else 1.0
            _enrich_rows.append({
                "phylon": _p,
                "kegg_pathway": _pw,
                "count": _k,
                "total_in_pathway": _K,
                "genes_in_phylon": _n,
                "p_value": _pval,
            })

    phylon_kegg_enrichment = pd.DataFrame(_enrich_rows)

    # Filter significant
    _sig = phylon_kegg_enrichment[phylon_kegg_enrichment["p_value"] < 0.05].copy()
    phylon_kegg_enrichment.to_csv(os.path.join(OUT, "data", "5c_phylon_kegg_enrichment.csv"), index=False)

    if _sig.empty:
        mo.output.replace(mo.md(f"No significant KEGG pathway enrichments (p < 0.05) found.\n\nSaved all {len(phylon_kegg_enrichment)} tests to CSV."))
        return (phylon_kegg_enrichment,)

    # Heatmap of significant pathways
    _sig_matrix = _sig.pivot_table(index="phylon", columns="kegg_pathway", values="p_value")
    _fig_kegg = sns.clustermap(_sig_matrix.fillna(0.05), cmap="rocket_r", figsize=(12, 6))
    _fig_kegg.savefig(os.path.join(FIG, "5c_kegg_heatmap.png"), bbox_inches="tight")
    plt.close(_fig_kegg.fig)

    mo.output.replace(
        mo.vstack([
            mo.md(
                f"KEGG pathway enrichments (p < 0.05):\n\n"
                f"- **{len(_sig)}** significant enrichments\n"
                f"- **{_sig['phylon'].nunique()}** phylons\n"
                f"- **{_sig['kegg_pathway'].nunique()}** pathways"
            ),
            mo.md("Saved `5c_kegg_heatmap.png`"),
        ])
    )
    return (phylon_kegg_enrichment,)


@app.cell
def _(L_BIN, phylon_cog_enrichment, phylon_go_enrichments, phylon_kegg_enrichment):
    """Summary of 5c outputs."""
    _lines = [
        "## 5c: Functional Enrichments — Summary\n",
        f"- L_binarized: {L_BIN.shape[0]} genes x {L_BIN.shape[1]} phylons\n",
        "**GO enrichments:**\n",
        f"- {len(phylon_go_enrichments)} significant (p < 0.05)\n",
        f"- {phylon_go_enrichments['phylon'].nunique()} phylons, "
        f"{phylon_go_enrichments['function'].nunique()} GO terms\n",
    ]
    if len(phylon_cog_enrichment) > 0:
        _sig_cog = phylon_cog_enrichment[phylon_cog_enrichment["p_value"] < 0.05]
        _lines.append("\n**COG enrichments:**\n")
        _lines.append(f"- {len(_sig_cog)} significant (p < 0.05)\n")
    if len(phylon_kegg_enrichment) > 0:
        _sig_kegg = phylon_kegg_enrichment[phylon_kegg_enrichment["p_value"] < 0.05]
        _lines.append("\n**KEGG enrichments:**\n")
        _lines.append(f"- {len(_sig_kegg)} significant (p < 0.05)\n")

    _lines.append(
        "\n**Outputs:**\n"
        "- `5c_phylon_go_enrichment.csv`, `5c_enrichment_heatmap.png`\n"
        "- `5c_phylon_cog_enrichment.csv`, `5c_cog_heatmap.png`\n"
        "- `5c_phylon_kegg_enrichment.csv`, `5c_kegg_heatmap.png`\n"
        "- `5c_wordcloud_phylon*.png` (one per phylon)"
    )
    mo.output.replace(mo.md("\n".join(_lines)))
    return


if __name__ == "__main__":
    app.run()
