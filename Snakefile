# Master Snakefile for pyphylon-marimo pipeline
# Orchestrates marimo notebook execution, bioinformatics workflows, and Quarto reports
#
# Usage:
#   snakemake -n                     # dry-run
#   snakemake --cores 4              # full pipeline
#   snakemake report_1a --cores 1    # single report

configfile: "config.yml"

TEMP = config["TEMP_DIR"]
OUT = config["OUTPUT_DIR"]
FIG = config["FIGURES_DIR"]
REP = config["REPORTS_DIR"]
DATA = config["DATA_DIR"]
SPECIES = config["PG_NAME"]

# =============================================================================
# Default target: all reports
# =============================================================================

rule all:
    input:
        # Phase 1 reports
        f"{REP}/1a_genome_filtering.pdf",
        f"{REP}/1b_genome_download.pdf",
        # Phase 2 reports
        f"{REP}/2a_clean_metadata.pdf",
        f"{REP}/2b_mash_filtration.pdf",
        f"{REP}/2c_build_pangenome.pdf",
        f"{REP}/2d_enrich_metadata.pdf",
        # Phase 3 reports
        f"{REP}/3a_extract_CAR.pdf",
        f"{REP}/3b_heaps_plot.pdf",
        # Phase 4 reports
        f"{REP}/4a_nmf_decomposition.pdf",
        # Phase 5 reports
        f"{REP}/5a_phylon_characterization.pdf",
        f"{REP}/5b_gene_diff.pdf",
        f"{REP}/5c_functional_enrichments.pdf",
        f"{REP}/5d_gene_alignment.pdf",
        f"{REP}/5e_blast_enrichment.pdf",
        f"{REP}/5f_infer_affinities.pdf",

# =============================================================================
# Phase 1: Genome Filtering & Download
# =============================================================================

rule nb_1a:
    input:
        config="config.yml"
    output:
        temp(f"{TEMP}/1a_genome_summary.csv"),
        temp(f"{TEMP}/1a_genome_metadata.csv"),
        f"{FIG}/1a_unfiltered_strains.png",
        f"{FIG}/1a_unfiltered_n50.png",
        f"{FIG}/1a_filtered_strains.png",
        f"{OUT}/data/1a_df_filtration.csv",
    shell:
        "python notebooks/1a_filter_genomes.py -- --config {input.config}"

rule nb_1b:
    input:
        config="config.yml",
        summary=f"{TEMP}/1a_genome_summary.csv",
        metadata=f"{TEMP}/1a_genome_metadata.csv",
    output:
        temp(f"{TEMP}/1b_genome_summary.csv"),
        temp(f"{TEMP}/1b_genome_metadata.csv"),
    shell:
        "python notebooks/1b_download_genomes.py -- --config {input.config}"

# =============================================================================
# Bioinformatics: Mash Workflow
# =============================================================================

# TODO: Integrate workflow/mash/Snakefile rules

# =============================================================================
# Phase 2: Mash Filtration & Pangenome Construction
# =============================================================================

rule nb_2a:
    input:
        config="config.yml",
        summary=f"{TEMP}/1b_genome_summary.csv",
        metadata=f"{TEMP}/1b_genome_metadata.csv",
    output:
        temp(f"{TEMP}/2a_genome_summary.csv"),
        temp(f"{TEMP}/2a_genome_metadata.csv"),
    shell:
        "python notebooks/2a_clean_metadata.py -- --config {input.config}"

rule nb_2b:
    input:
        config="config.yml",
        summary=f"{TEMP}/2a_genome_summary.csv",
        metadata=f"{TEMP}/2a_genome_metadata.csv",
        mash_dist=f"{TEMP}/2b_mash/mash_distances.txt",
    output:
        temp(f"{TEMP}/2b_genome_summary.csv"),
        temp(f"{TEMP}/2b_genome_metadata.csv"),
        f"{OUT}/data/2b_mash_square.csv",
        f"{OUT}/data/2b_mash_corr_dist.csv",
        f"{FIG}/2b_mash_heatmap.png",
        f"{FIG}/2b_mash_dist_all.png",
        f"{FIG}/2b_mash_dist_filtered.png",
        f"{FIG}/2b_sensitivity_analysis.png",
        f"{FIG}/2b_clustermap_initial.png",
        f"{FIG}/2b_clustermap_final.png",
        f"{FIG}/2b_cluster_sizes_initial.png",
        f"{FIG}/2b_cluster_sizes_final.png",
    shell:
        "python notebooks/2b_mash_filtration.py -- --config {input.config}"

rule nb_2c:
    input:
        config="config.yml",
        metadata=f"{TEMP}/2b_genome_metadata.csv",
    output:
        f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
        f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_allele.pickle.gz",
        f"{FIG}/2c_genes_per_genome.png",
        f"{OUT}/data/2c_pangenome_summary.csv",
    shell:
        "python notebooks/2c_build_pangenome.py -- --config {input.config}"

rule nb_2d:
    input:
        config="config.yml",
        metadata=f"{TEMP}/2b_genome_metadata.csv",
        mlst=f"{DATA}/processed/mlst_report.txt",
    output:
        temp(f"{TEMP}/2d_enriched_metadata.csv"),
    shell:
        "python notebooks/2d_enrich_metadata.py -- --config {input.config}"

# =============================================================================
# Bioinformatics: Annotation, MLST, CD-HIT Workflow
# =============================================================================

# TODO: Integrate workflow/anno_mlst_cdhit/Snakefile rules

# =============================================================================
# Phase 3: Genome Characterization
# =============================================================================

rule nb_3a:
    input:
        config="config.yml",
        pickle=f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
        metadata=f"{TEMP}/2d_enriched_metadata.csv",
    output:
        f"{DATA}/processed/CAR_genomes/df_core.csv",
        f"{DATA}/processed/CAR_genomes/df_acc.csv",
        f"{DATA}/processed/CAR_genomes/df_rare.csv",
        f"{FIG}/3a_gene_frequency.png",
        f"{FIG}/3a_pangenome_segments.png",
        f"{OUT}/data/3a_car_summary.csv",
    shell:
        "python notebooks/3a_extract_CAR.py -- --config {input.config}"

rule nb_3b:
    input:
        config="config.yml",
        pickle=f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
        metadata=f"{TEMP}/2d_enriched_metadata.csv",
    output:
        f"{FIG}/3b_heaps_law.png",
        f"{OUT}/data/3b_heaps_summary.csv",
    shell:
        "python notebooks/3b_heaps_plot.py -- --config {input.config}"

# =============================================================================
# Phase 4: NMF Decomposition
# =============================================================================

# TODO: Add rule for 4a

# =============================================================================
# Phase 5: Phylon Analysis
# =============================================================================

# TODO: Add rules for 5a-5f

# =============================================================================
# Bioinformatics: Infer Affinities Workflow
# =============================================================================

# TODO: Integrate workflow/infer_affinities/Snakefile rules

# =============================================================================
# Quarto Reports
# =============================================================================

rule report_1a:
    input:
        f"{FIG}/1a_unfiltered_strains.png",
        f"{FIG}/1a_unfiltered_n50.png",
        f"{FIG}/1a_filtered_strains.png",
        f"{OUT}/data/1a_df_filtration.csv",
    output:
        f"{REP}/1a_genome_filtering.pdf"
    shell:
        "quarto render reports/1a_genome_filtering.qmd --to pdf"

rule report_1b:
    input:
        f"{TEMP}/1b_genome_summary.csv",
        f"{TEMP}/1b_genome_metadata.csv",
    output:
        f"{REP}/1b_genome_download.pdf"
    shell:
        "quarto render reports/1b_genome_download.qmd --to pdf"

rule report_2a:
    input:
        f"{TEMP}/2a_genome_summary.csv",
        f"{TEMP}/2a_genome_metadata.csv",
    output:
        f"{REP}/2a_clean_metadata.pdf"
    shell:
        "quarto render reports/2a_clean_metadata.qmd --to pdf"

rule report_2b:
    input:
        f"{TEMP}/2b_genome_metadata.csv",
        f"{OUT}/data/2b_mash_square.csv",
        f"{OUT}/data/2b_mash_corr_dist.csv",
        f"{FIG}/2b_mash_heatmap.png",
        f"{FIG}/2b_mash_dist_all.png",
        f"{FIG}/2b_mash_dist_filtered.png",
        f"{FIG}/2b_sensitivity_analysis.png",
        f"{FIG}/2b_clustermap_initial.png",
        f"{FIG}/2b_clustermap_final.png",
        f"{FIG}/2b_cluster_sizes_initial.png",
        f"{FIG}/2b_cluster_sizes_final.png",
    output:
        f"{REP}/2b_mash_filtration.pdf"
    shell:
        "quarto render reports/2b_mash_filtration.qmd --to pdf"

rule report_2c:
    input:
        f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
        f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_allele.pickle.gz",
        f"{FIG}/2c_genes_per_genome.png",
        f"{OUT}/data/2c_pangenome_summary.csv",
    output:
        f"{REP}/2c_build_pangenome.pdf"
    shell:
        "quarto render reports/2c_build_pangenome.qmd --to pdf"

rule report_2d:
    input:
        f"{TEMP}/2d_enriched_metadata.csv",
    output:
        f"{REP}/2d_enrich_metadata.pdf"
    shell:
        "quarto render reports/2d_enrich_metadata.qmd --to pdf"

rule report_3a:
    input:
        f"{DATA}/processed/CAR_genomes/df_core.csv",
        f"{DATA}/processed/CAR_genomes/df_acc.csv",
        f"{DATA}/processed/CAR_genomes/df_rare.csv",
        f"{FIG}/3a_gene_frequency.png",
        f"{FIG}/3a_pangenome_segments.png",
        f"{OUT}/data/3a_car_summary.csv",
    output:
        f"{REP}/3a_extract_CAR.pdf"
    shell:
        "quarto render reports/3a_extract_CAR.qmd --to pdf"

rule report_3b:
    input:
        f"{FIG}/3b_heaps_law.png",
        f"{OUT}/data/3b_heaps_summary.csv",
    output:
        f"{REP}/3b_heaps_plot.pdf"
    shell:
        "quarto render reports/3b_heaps_plot.qmd --to pdf"

# TODO: Add report rules as notebooks are converted
