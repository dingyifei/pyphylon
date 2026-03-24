import marimo

__generated_with = "0.20.4"
app = marimo.App(
    width="medium",
    app_title="Presentation",
    layout_file="layouts/presentation.slides.json",
)

with app.setup:
    import os
    import re
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

    # Load phylon characterization data
    df_profiles = pd.read_csv(os.path.join(DATA, "5a_phylon_profiles.csv"))
    df_strain_demo = pd.read_csv(os.path.join(DATA, "5a_phylon_strain_demographics.csv"))

    # Load VFDB cross-reference (if available)
    _vfdb_path = os.path.join(DATA, "vfdb_phylon_category_summary.csv")
    df_vfdb = pd.read_csv(_vfdb_path) if os.path.exists(_vfdb_path) else None
    return FIG, SPECIES, df_car, df_filt, df_vfdb, nmf, pg, phy


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

    4 mobile genetic elements &bull; 9 lineage-defining clonal groups
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
            mo.image(src=os.path.join(FIG, "1a_unfiltered_strains.png"), width="50%"),
        ]
    )
    return


@app.cell
def _(FIG, df_filt):
    _after_n50 = int(df_filt.loc["L50/N50", "remaining"])
    _after_checkm = int(df_filt.loc["CheckM_completeness_contamination", "remaining"])
    mo.vstack(
        [
            mo.md(
                f"""
        ## Filtering

        | Filter | Genomes remaining |
        |--------|------------------:|
        | L50 / N50 threshold | {_after_n50} |
        | CheckM completeness &ge; 92%, contamination &lt; 5% | {_after_checkm} |
        | Mash species-level verification | **465** |
        """
            ),
            mo.hstack(
                [
                    mo.image(src=os.path.join(FIG, "1a_unfiltered_n50.png"), width="50%"),
                    mo.image(src=os.path.join(FIG, "1a_filtered_strains.png"), width="50%"),
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
            mo.image(src=os.path.join(FIG, "2b_mash_heatmap.png"), width="50%"),
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
            mo.image(src=os.path.join(FIG, "2b_clustermap_final.png"), width="50%"),
        ]
    )
    return


@app.cell
def _():
    mo.md("""
    ## Annotation Engine Update: Bakta 1.12.0 + DB v6.0

    Problem: lots of hypothetical proteins

    TLDR:

    DB v6.0 brings substantially updated reference databases:

    | Database | v5.1 | v6.0 | Impact |
    |----------|------|------|--------|
    | **UniProtKB** | 2023_05 | **2025_01** | Fewer hypothetical proteins, better GO terms |
    | **RefSeq** | r221 | **r228** | 7 releases of new prokaryotic annotations |
    | **AMRFinderPlus** | 2023-11 | **2024-12** | More complete AMR gene detection |
    | **COG** | 2020 | **2024** | 4 years of updated functional categories |
    | **Pfam** | 36 | **37.2** | New domain families, improved HMM profiles |
    | **VFDB** | 2024-01 | **2025-02** | Expanded virulence factor catalog |
    | **Rfam** | 14.10 | **15** | Updated ncRNA families |
    """)
    return


@app.cell
def _():
    """Annotation coverage improvement from eggNOG-mapper + COGclassifier."""
    _total = 8096  # representative sequences

    # Load eggNOG annotations
    _eg = pd.read_csv(
        os.path.join("data", "processed", "eggnog", "CJejuni.emapper.annotations"),
        sep="\t", comment="#", header=None,
    )
    _eg_clusters = _eg[0].apply(lambda x: re.sub(r"A\d+$", "", x))
    _eg_cog = (_eg[6] != "-").sum()
    _eg_kegg = (_eg[11] != "-").sum()

    # Load COGclassifier
    _cog_cl = pd.read_csv(
        os.path.join("data", "processed", "cogclassifier", "cog_classify.tsv"), sep="\t",
    )
    _cog_cl_clusters = set(_cog_cl["QUERY_ID"].apply(lambda x: re.sub(r"A\d+$", "", x)))
    _eg_cog_clusters = set(_eg_clusters[_eg[6] != "-"])
    _union_cog = len(_eg_cog_clusters | _cog_cl_clusters)

    mo.md(f"""
    ## Annotation Coverage: Beyond BAKTA

    BAKTA db-full provides product names but limited **functional classification**.
    We added **eggNOG-mapper** (orthology-based) and **COGclassifier** (NCBI RPS-BLAST):

    | Annotation | BAKTA db-full | + eggNOG-mapper | + COGclassifier | Combined |
    | ---------- | ------------- | --------------- | --------------- | -------- |
    | **COG categories** | <1% | {_eg_cog * 100 / _total:.1f}% ({_eg_cog:,}) | {len(_cog_cl_clusters) * 100 / _total:.1f}% ({len(_cog_cl_clusters):,}) | **{_union_cog * 100 / _total:.1f}%** ({_union_cog:,}) |
    | **KEGG KO** | ~29% | {_eg_kegg * 100 / _total:.1f}% ({_eg_kegg:,}) | n/a | ~55%+ |

    *{_total:,} representative protein sequences (one per CD-HIT cluster)*

    **COG coverage jumped from <1% to {_union_cog * 100 / _total:.0f}%** — enabling
    functional category distribution analysis per phylon (comparable to *E. coli* study).
    """)
    return


@app.cell
def _(FIG, pg):
    mo.vstack(
        [
            mo.md(
                f"""
        ## Pangenome Construction

        Genomes annotated with **Bakta 1.12.0**, protein sequences clustered with
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

        power-law model fitting:

        | Category | Genes | % of pangenome |
        |----------|------:|---------------:|
        | **Core** (present in &gt;98.5% of strains) | {_core:,} | {100 * _core / _total:.1f}% |
        | **Accessory** | {_acc:,} | {100 * _acc / _total:.1f}% |
        | **Rare** (present in &lt;8.8% of strains) | {_rare:,} | {100 * _rare / _total:.1f}% |
        | **Total** | {_total:,} | 100% |
        """
            ),
            mo.image(src=os.path.join(FIG, "3a_pangenome_segments.png"), width="50%"),
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
            mo.image(src=os.path.join(FIG, "3a_gene_frequency.png"), width="50%"),
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

        - Best rank: **{int(nmf["best_rank"])}**
        - Best AIC: **{nmf["best_aic"]:,.1f}**
        """
            ),
            mo.image(src=os.path.join(FIG, "4a_mca_variance.png"), width="50%"),
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
            mo.image(src=os.path.join(FIG, "4a_consensus_clustermap.png"), width="50%"),
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
            mo.image(src=os.path.join(FIG, "5a_regression.png"), width="50%"),
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
            mo.vstack(
                [
                    mo.image(src=os.path.join(FIG, "5a_L_binarized_sorted.png"), width="50%"),
                    mo.image(src=os.path.join(FIG, "5a_A_binarized_sorted.png"), width="50%"),
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
        ## Phylon Dendrogram: MGE vs Lineage Split

        Hierarchical clustering reveals a **bipartite architecture**:

        - **Root split** perfectly separates **4 MGE phylons** {4, 5, 7, 10}
          from **9 lineage-defining phylons** {0, 1, 2, 3, 6, 8, 9, 11, 12}
        - MGE phylons = phages + conjugative elements (horizontally transferred)
        - Lineage phylons = clonal population structure (vertically inherited)
        - **Split 2**: Phylon 1 (Finnish ST-677) separates with 93 exclusive genes
        - **Split 7**: CC-21-related cluster {0, 2, 8} groups together
        """
            ),
            mo.image(src=os.path.join(FIG, "5b_phylon_dendrogram.png"), width="50%"),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Four Mobile Genetic Element Phylons

        | Phylon | Element | Genes | Strains | Key features |
        |--------|---------|------:|--------:|-------------|
        | **4** | Mu-like prophage | 51 | 130 (28%) | Transposable phage, P2 tail, integrase, holin/lysozyme |
        | **5** | T4SS conjugative | 40 | 114 (25%) | pVir/pTet-like, VirB4/B9/B10/B11, **Tet(O) resistance** |
        | **7** | CJIE1-like phage | 47 | 52 (11%) | HK97 temperate siphovirus, Gam nuclease inhibitor |
        | **10** | Lambda-like phage | 34 | 65 (14%) | Bet/YqaJ recombination, RusA endonuclease |

        - Phylon 5 has **100% exclusive genes** -- every gene is unique to this element
        - Phylon 4 (Mu-like) is the **most widespread** MGE across all lineages
        - Phylon 5 carries **tetracycline resistance** via Tet(O) ribosomal protection
        - Conjugation is **thermoregulated at 42&deg;C** (avian body temperature)
        """
            ),
            mo.image(src=os.path.join(FIG, "5b_first_split_heatmap.png")),
        ]
    )
    return


@app.cell
def _():
    mo.md("""
    ## Nine Lineage-Defining Phylons

    | Phylon | MLST | Geography | Host | Signature |
    |--------|------|-----------|------|-----------|
    | **0** | ST-11422 | USA (91%) | Unknown | R-M system, transposases |
    | **1** | ST-677 | Finland (100%) | Human (100%) | 158 exclusive genes, phage remnants |
    | **2** | ST-4253/CC-21 | S. Korea (57%) | Mixed (incl. dog) | CPS biosynthesis |
    | **3** | ST-45/CC-45 | Cosmopolitan | Mixed | Type I R-M, lactate permease |
    | **6** | ST-2993 | Peru (90%) | Human (100%) | Glycosyltransferases, **GBS outbreak** |
    | **8** | ST-50/CC-21 | USA (40%) | Chicken (24%) | Sialyltransferase, LOS biosynthesis |
    | **9** | Diverse | USA/India | Mixed | **T6SS**, PEB1 adhesin |
    | **11** | ST-1044 | S. Korea (39%) | Chicken (45%) | LOS glycosyltransferase |
    | **12** | ST-22/CC-22 | NL/USA | Human (68%) | **HS:19, Guillain-Barr&eacute; syndrome** |
    """)
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Phylon 12: HS:19 / Guillain-Barr&eacute; Syndrome Lineage

        - **ST-22 (72%)**, serovar **HS:19 (60%)** -- strongest GBS association
        - LOS ganglioside mimicry (GM1/GD1a/GQ1b) triggers autoimmune neuropathy
        - Sialylated LOS class A1 with *cst-II* sialyltransferase
        - Netherlands + USA distribution (Atlantic corridor)
        - Phylon 6 (ST-2993, Peru) = **second independent GBS lineage**
          linked to 2019 Peru outbreak
        """
            ),
            mo.image(src=os.path.join(FIG, "5d_circular_phylon12.png")),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## Phylon 1: Finnish ST-677 Bacteremia Clone

        - **100% Finland**, 100% human, 100% mash cluster 9
        - **Most exclusive genes** (158) of any lineage phylon
        - Associated with **invasive bacteremia** (50% of Finnish blood isolates)
        - Serum resistance via unique CPS locus (*C. doylei*-like)
        - Lost 93 conserved genes present in all other lineages
        - CDT operon degradation -- unique virulence profile
        """
            ),
            mo.image(src=os.path.join(FIG, "5d_circular_phylon1.png")),
        ]
    )
    return


@app.cell
def _(FIG):
    mo.vstack(
        [
            mo.md(
                """
        ## GO Functional Enrichment

        Hypergeometric testing of phylon gene sets against GO annotations
        reveals functional specialization:

        - **CPS/LOS biosynthesis** enriched in most lineage phylons
        - **DNA integration/recombination** enriched in MGE phylons
        - Phylon-specific functions: lactate transport (phylon 3),
          sialic acid biosynthesis (phylon 8), secretion (phylon 9)
        """
            ),
            mo.image(src=os.path.join(FIG, "5c_enrichment_heatmap.png")),
        ]
    )
    return


@app.cell
def _(df_vfdb):
    if df_vfdb is not None:
        mo.md("""
        ## VFDB Virulence Factor Cross-Reference

        434 BLAST hits against VFDB mapped to phylon gene assignments:

        | Category | Hits | Key finding |
        |----------|-----:|------------|
        | Motility (flagella) | 212 | Distributed across all lineage phylons |
        | Immune modulation (CPS/LOS) | 187 | Highest in CC-21 phylons (2, 8) |
        | Exotoxin (CDT) | 16 | Present in 7 of 9 lineage phylons |
        | Adherence (PEB1) | 10 | **Exclusive to phylon 9** |
        | Effector delivery | 6 | T4SS-associated (phylon 5) |
        | Invasion (CiaB) | 3 | Phylon 5 carries CiaB on conjugative element |

        - MGE phylons have **0-3 VFDB hits** (mobile elements, not virulence loci)
        - 371 of 434 hits are **core genes** (not in accessory set)
        """)
    else:
        mo.md("## VFDB Cross-Reference\n\n*VFDB summary not yet generated.*")
    return


@app.cell
def _():
    mo.md("""
    ## Geographic Structuring of Phylons

    | Pattern | Phylons | Interpretation |
    |---------|---------|---------------|
    | **Finland-exclusive** | 1 (ST-677) | Locally adapted bacteremia clone |
    | **Peru-restricted** | 6 (ST-2993) | Geographically isolated GBS lineage |
    | **East Asian** | 2, 11 (S. Korea + China) | Regional poultry lineages |
    | **North American** | 0, 7 (USA/Canada) | USA-dominant |
    | **Cosmopolitan** | 3, 8, 9 (CC-45, CC-21) | Globally distributed generalists |
    | **Atlantic** | 12 (NL + USA) | GBS-associated corridor |

    Strong geographic signals reflect both genuine population structure
    and sampling bias in the BV-BRC database.
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Host Adaptation Signatures

    **Poultry-associated** (&gt;20% chicken):
    - Phylon 8 (CC-21): Sialylation machinery, aerotolerance, cold tolerance
    - Phylon 11 (East Asian): LOS glycosyltransferase, motility factors

    **Human-exclusive** (no known animal reservoir):
    - Phylon 1 (ST-677): Finnish bacteremia clone, serum resistance
    - Phylon 6 (ST-2993): Peruvian GBS lineage

    **Generalist:**
    - Phylon 3 (CC-45): Humans, chickens, bears, dogs -- rapid host switching
    - Phylon 9 (diverse): T6SS + PEB1 adhesin across humans, chickens, cattle

    CC-21 (phylon 8) and CC-45 (phylon 3) are both globally dominant generalists
    but use **divergent strategies**: CC-21 is stress-adapted with AMR accumulation,
    CC-45 is pansusceptible with distinct metabolism.
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Key Findings

    1. **Binary MGE-lineage architecture** -- root split perfectly separates
       4 mobile elements from 9 clonal lineages, validating NMF decomposition

    2. **Two independent GBS lineages** -- phylon 12 (HS:19/ST-22, classic)
       and phylon 6 (ST-2993, 2019 Peru outbreak) evolved neuropathogenicity independently

    3. **Phylon 5 = tetracycline resistance vehicle** -- pVir/pTet conjugative
       element in 24.5% of strains, thermoregulated transfer at avian 42&deg;C

    4. **CC-677 bacteremia clone** (phylon 1) -- 158 exclusive genes, serum
       resistance, unique CPS, lost 93 conserved genes

    5. **CPS diversity drives lineage differentiation** -- immune modulation
       (CPS/LOS) is the top functional enrichment for 8 of 13 phylons

    6. **T6SS discovery** in phylon 9 -- unusual for *C. jejuni*, with
       PEB1 adhesin exclusivity suggesting novel adhesion strategies
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Conclusions

    1. ***C. jejuni*'s open pangenome** (Heaps' &gamma; &asymp; 0.28) decomposes
       cleanly into **13 biologically interpretable phylons** (R&sup2; = 0.94)

    2. **NMF captures both vertical and horizontal inheritance** -- lineage
       phylons map to known clonal complexes; MGE phylons map to characterized
       phages and conjugative elements

    3. **Clinically actionable insights**: GBS risk lineages (phylons 6, 12),
       antibiotic resistance vehicle (phylon 5), bacteremia clone (phylon 1)

    4. **First NMF phylon decomposition of *C. jejuni*** -- method previously
       applied only to *E. coli* and *Enterobacter* (Palsson lab)

    5. **Fully reproducible** pipeline: Snakemake + marimo notebooks, from
       raw genome download through biological characterization
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Methods Summary

    | Stage | Tool / Method |
    |-------|---------------|
    | Genome retrieval | BV-BRC API (complete genomes only) |
    | Quality control | CheckM (completeness &ge; 92%, contamination &lt; 5%) |
    | Species verification | Mash distance clustering |
    | Genome annotation | Bakta 1.12.0 + DB v6.0 |
    | Sequence typing | MLST |
    | Pangenome construction | CD-HIT (80% identity, 80% coverage) |
    | Gene classification | Power-law model (Core/Accessory/Rare) |
    | Matrix factorization | NMF with AIC rank selection (k=13) |
    | Phylon characterization | GO enrichment, VFDB cross-ref, circular genome plots |
    | Literature context | PubMed systematic search (34 references) |
    | Pipeline orchestration | Snakemake + marimo notebooks |
    """)
    return


if __name__ == "__main__":
    app.run()
