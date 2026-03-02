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

    from pyphylon.biointerp import (
        calc_all_phylon_go_enrichments,
        collect_functions,
        explode_go_annos,
        gen_phylon_wordcloud,
        get_go_mapping,
        get_pg_to_locus_map,
    )

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
    SPECIES = CONFIG["PG_NAME"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)

    return CONFIG, DATA, FIG, OUT, SPECIES


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
def _(FIG, L_BIN, OUT, cluster_to_go_functions, functions2genes, go_functions):
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
def _(FIG, phylon_go_enrichments):
    """Generate enrichment heatmap (phylon × GO term, colored by p-value)."""
    phylon_go_enrichments_mat = pd.pivot_table(
        phylon_go_enrichments, index="phylon", columns="function", values="p_value"
    )
    g = sns.clustermap(phylon_go_enrichments_mat.fillna(0.05), cmap="rocket_r")
    g.savefig(os.path.join(FIG, "5c_enrichment_heatmap.png"), bbox_inches="tight")
    plt.close(g.fig)

    mo.output.replace(mo.md("Saved `5c_enrichment_heatmap.png`"))


@app.cell
def _(FIG, L_BIN, functions2genes, phylon_go_enrichments):
    """Generate and save wordclouds for each phylon."""
    wordcloud_phylons = sorted(phylon_go_enrichments["phylon"].unique())
    for _phylon in wordcloud_phylons:
        wc_path = os.path.join(FIG, f"5c_wordcloud_{_phylon}.png")
        gen_phylon_wordcloud(L_BIN, functions2genes, _phylon, cutoff=0, save=True, filename=wc_path)

    mo.output.replace(mo.md(f"Saved **{len(wordcloud_phylons)}** wordcloud PNGs to `output/figures/`"))


@app.cell
def _(L_BIN, phylon_go_enrichments):
    """Summary of 5c outputs."""
    mo.output.replace(
        mo.md(
            f"## 5c: Functional Enrichments — Summary\n\n"
            f"- L_binarized: {L_BIN.shape[0]} genes x {L_BIN.shape[1]} phylons\n"
            f"- Significant enrichments: {len(phylon_go_enrichments)} (p < 0.05)\n"
            f"- Phylons with enrichments: {phylon_go_enrichments['phylon'].nunique()}\n"
            f"- Unique GO terms: {phylon_go_enrichments['function'].nunique()}\n\n"
            f"**Outputs:**\n"
            f"- `5c_phylon_go_enrichment.csv`\n"
            f"- `5c_enrichment_heatmap.png`\n"
            f"- `5c_wordcloud_phylon*.png` (one per phylon)"
        )
    )


if __name__ == "__main__":
    app.run()
