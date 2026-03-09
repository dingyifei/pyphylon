import marimo

__generated_with = "0.20.2"
app = marimo.App(
    width="medium",
    app_title="Presentation",
    layout_file="layouts/presentation.slides.json",
)

with app.setup:
    import os
    import sys

    import marimo as mo
    import pandas as pd
    import yaml


@app.cell
def _():
    """Parse config and load all summary data."""

    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    FIG = CONFIG["FIGURES_DIR"]
    DATA = os.path.join(CONFIG["OUTPUT_DIR"], "data")
    SPECIES = CONFIG["SPECIES_NAME"]

    # Load summary CSVs
    df_filt = pd.read_csv(os.path.join(DATA, "1a_df_filtration.csv"), index_col=0)
    df_pg = pd.read_csv(os.path.join(DATA, "2c_pangenome_summary.csv"))
    pg = dict(zip(df_pg["metric"], df_pg["value"]))

    df_car = pd.read_csv(os.path.join(DATA, "3a_car_summary.csv"))

    df_nmf = pd.read_csv(os.path.join(DATA, "4a_nmf_summary.csv"))
    nmf = dict(zip(df_nmf["metric"], df_nmf["value"]))

    df_phy = pd.read_csv(os.path.join(DATA, "5a_phylon_summary.csv"))
    phy = dict(zip(df_phy["metric"], df_phy["value"]))
    return FIG, SPECIES, df_car, df_filt, nmf, pg, phy


@app.cell
def _():
    mo.md("""
    ## Background: What is *Campylobacter jejuni*?

    - **Leading cause of bacterial gastroenteritis** worldwide,
      responsible for an estimated 550 million cases per year (WHO)
    - Gram-negative, microaerophilic, spiral-shaped bacterium
    - Commonly found in poultry and other animal reservoirs
    - **Highly recombinogenic** — frequent HGT
    - Infections can trigger Guillain-Barre syndrome, a serious autoimmune neuropathy
    - Phylons may reveal functional modules linked to host adaptation, virulence,
      and antimicrobial resistance
    """)
    return


@app.cell
def _(SPECIES, pg, phy):
    mo.md(f"""
    # Pangenomic Analysis of *{SPECIES}*

    ## Identifying Phylons via Non-negative Matrix Factorization

    ---

    **{int(pg["n_genomes"])}** complete genomes &bull;
    **{int(pg["n_genes"]):,}** gene families &bull;
    **{int(phy["n_phylons"])}** phylons identified

    TODO: ADD filter stats
    """)
    return


@app.cell
def _(FIG, df_filt):
    _initial = int(df_filt.loc["prefiltration", "remaining"])
    mo.vstack(
        [
            mo.md(
                f"""
        ## Genome Selection

        Queried the [BV-BRC](https://www.bv-brc.org/) database for **complete**,
        **good-quality** *C. jejuni* genomes.

        **{_initial} genomes** retrieved — genome length vs. gene count shown below.
        """
            ),
            mo.image(src=os.path.join(FIG, "1a_unfiltered_strains.png")),
        ]
    )

    # TODO: switch to the plot after filter
    return


@app.cell
def _(FIG, df_filt):
    _after_n50 = int(df_filt.loc["L50/N50", "remaining"])
    _after_checkm = int(df_filt.loc["CheckM_completeness_contamination", "remaining"])
    mo.vstack(
        [
            mo.md(
                f"""
        ## Quality Control Filtering

        | Filter | Genomes remaining |
        |--------|------------------:|
        | L50 / N50 threshold | {_after_n50} |
        | CheckM completeness &ge; 92%, contamination &lt; 5% | {_after_checkm} |
        | Mash species-level verification | **465** |
        """
            ),
            mo.hstack(
                [
                    mo.image(src=os.path.join(FIG, "1a_unfiltered_n50.png")),
                    mo.image(src=os.path.join(FIG, "1a_filtered_strains.png")),
                ]
            ),
        ]
    )
    # TODO: add checkM plot, mention the filter table
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Mash Distance Analysis

        Pairwise genome distances computed via **Mash** (MinHash sketching).
        Hierarchical clustering with Ward's linkage identifies strain groups
        and removes outliers.
        """
            ),
            mo.image(src=os.path.join(FIG, "2b_mash_heatmap.png")),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Mash Clustering Result

        Ward's linkage dendrogram with cluster annotations.
        Elbow-based threshold selection balances cluster granularity
        and biological interpretability.
        """
            ),
            mo.image(src=os.path.join(FIG, "2b_clustermap_final.png")),
        ]
    )
    return


@app.cell
def _(FIG, pg):
    mo.vstack(
        [
            mo.md(
                f"""
        ## Pangenome Construction

        Genomes annotated with **Bakta**, protein sequences clustered with
        **CD-HIT** (80% identity, 80% alignment coverage).

        | Metric | Value |
        |--------|------:|
        | Genomes | {int(pg["n_genomes"])} |
        | Gene families | {int(pg["n_genes"]):,} |
        | Total alleles | {int(pg["n_alleles"]):,} |
        | Mean genes/genome | {pg["mean_genes_per_genome"]:.1f} &pm; {pg["std_genes_per_genome"]:.1f} |
        | Range | {int(pg["min_genes_per_genome"])} &ndash; {int(pg["max_genes_per_genome"])} |
        """
            ),
            mo.image(src=os.path.join(FIG, "2c_genes_per_genome.png")),
        ]
    )
    return


@app.cell
def _(FIG, df_car):
    _core = df_car[df_car["Category"] == "Core"]["Genes"].values[0]
    _acc = df_car[df_car["Category"] == "Accessory"]["Genes"].values[0]
    _rare = df_car[df_car["Category"] == "Rare"]["Genes"].values[0]
    _total = df_car[df_car["Category"] == "Total"]["Genes"].values[0]
    mo.vstack(
        [
            mo.md(
                f"""
        ## Core / Accessory / Rare Genes

        Gene families partitioned by frequency across strains using
        power-law model fitting:

        | Category | Genes | % of pangenome |
        |----------|------:|---------------:|
        | **Core** (present in &gt;98.5% of strains) | {_core:,} | {100 * _core / _total:.1f}% |
        | **Accessory** | {_acc:,} | {100 * _acc / _total:.1f}% |
        | **Rare** (present in &lt;8.8% of strains) | {_rare:,} | {100 * _rare / _total:.1f}% |
        | **Total** | {_total:,} | 100% |
        """
            ),
            mo.image(src=os.path.join(FIG, "3a_pangenome_segments.png")),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Open Pangenome

        **Heaps' law** models pangenome growth: G(n) = &kappa; &middot; n<sup>&gamma;</sup>

        - Coefficient &gamma; &asymp; **0.28** &lt; 1 &rarr; **open pangenome**
        - Each new genome continues to contribute novel genes
        - Consistent with *C. jejuni*'s high recombination rate and
          ecological diversity
        """
            ),
            mo.image(src=os.path.join(FIG, "3a_gene_frequency.png")),
        ]
    )
    return


@app.cell
def _(FIG, nmf):
    mo.vstack(
        [
            mo.md(
                f"""
        ## NMF Rank Selection

        **Multiple Correspondence Analysis (MCA)** provides an initial
        dimensionality estimate. NMF is then run across multiple ranks
        and the best is selected by **Akaike Information Criterion (AIC)**.

        - Best rank: **{int(nmf["best_rank"])}**
        - Best AIC: **{nmf["best_aic"]:,.1f}**
        """
            ),
            mo.image(src=os.path.join(FIG, "4a_mca_variance.png")),
        ]
    )

    # Missing the NMF decomp table
    return


@app.cell
def _(FIG, nmf):
    mo.vstack(
        [
            mo.md(
                f"""
        ## Consensus Clustering

        Multiple NMF runs at the best rank are aggregated into a
        **consensus matrix**. Clusters are identified via hierarchical
        clustering on the consensus.

        | Metric | Full | Filtered |
        |--------|-----:|---------:|
        | Clusters | {int(nmf["n_clusters_full"])} | {int(nmf["n_clusters_filtered"])} |
        | Cophenetic correlation | {nmf["cophenetic_full"]:.4f} | {nmf["cophenetic_filtered"]:.4f} |
        """
            ),
            mo.image(src=os.path.join(FIG, "4a_consensus_clustermap.png")),
        ]
    )
    return


@app.cell
def _(FIG, phy):
    mo.vstack(
        [
            mo.md(
                f"""
        ## Phylon Identification

        NMF decomposition refined into discrete phylon memberships:

        - **{int(phy["n_phylons"])} phylons** identified
        - **{int(phy["n_genes"]):,} genes** assigned to at least one phylon
        - **{int(phy["n_genomes"])} genomes** with phylon membership
        - Regression R&sup2; = **{phy["regression_r2"]:.3f}**
          (gene frequency vs. multi-phylon membership)
        """
            ),
            mo.image(src=os.path.join(FIG, "5a_regression.png")),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Binarized Phylon Matrices

        **L matrix** (strains &times; phylons) binarized via 3-means clustering.
        **A matrix** (phylons &times; genes) binarized at threshold &ge; 0.5.

        Block-diagonal structure confirms well-separated phylon assignments.
        """
            ),
            mo.hstack(
                [
                    mo.image(src=os.path.join(FIG, "5a_L_binarized_sorted.png")),
                    mo.image(src=os.path.join(FIG, "5a_A_binarized_sorted.png")),
                ]
            ),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Gene Differentiation

        Hierarchical clustering of phylons by shared gene content reveals
        the **phylon hierarchy** — groups of phylons that share large gene
        blocks versus those with distinct, exclusive gene repertoires.
        """
            ),
            mo.image(src=os.path.join(FIG, "5b_phylon_dendrogram.png")),
        ]
    )
    return


@app.cell
def _(phy):
    mo.md(f"""
    ## Conclusions

    1. ***C. jejuni* has an open pangenome** — new genomes continue to
       add novel genes (Heaps' &gamma; &asymp; 0.28), reflecting high
       recombination and ecological diversity

    2. **NMF successfully decomposed the pangenome into
       {int(phy["n_phylons"])} phylons** with a strong factorization
       fit (R&sup2; = {phy["regression_r2"]:.3f})

    3. **Phylons represent biologically meaningful modules** —
       co-occurring gene sets that may correspond to mobile elements,
       metabolic islands, or lineage-specific adaptations

    4. **The pipeline is fully reproducible** via Snakemake + marimo
       notebooks, from raw genome download to phylon characterization
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Current Methods Summary

    | Stage | Tool / Method |
    |-------|---------------|
    | Genome retrieval | BV-BRC API (complete genomes only) |
    | Quality control | CheckM (completeness &ge; 92%, contamination &lt; 5%) |
    | Species verification | Mash distance clustering |
    | Genome annotation | Bakta |
    | Sequence typing | MLST |
    | Pangenome construction | CD-HIT (80% identity, 80% coverage) |
    | Gene classification | Power-law model fitting (Core/Accessory/Rare) |
    | Dimensionality reduction | Multiple Correspondence Analysis (MCA) |
    | Matrix factorization | Non-negative Matrix Factorization (NMF) |
    | Consensus clustering | Hierarchical clustering on consensus matrix |
    | Binarization | 3-means clustering (L), threshold &ge; 0.5 (A) |
    | Pipeline orchestration | Snakemake + marimo notebooks |
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## This week's Tasks

    - []
    """)
    return


if __name__ == "__main__":
    app.run()
