# Master Snakefile for pyphylon-marimo pipeline
# Orchestrates marimo notebook execution and bioinformatics workflows
#
# Usage:
#   snakemake -n                     # dry-run
#   snakemake --cores 4              # full pipeline
#   snakemake nb_1a --cores 1        # single notebook

configfile: "config.yml"

TEMP = config["TEMP_DIR"]
OUT = config["OUTPUT_DIR"]
FIG = config["FIGURES_DIR"]
DATA = config["DATA_DIR"]
SPECIES = config["PG_NAME"]
FILTERED_GENOMES_FILE = config["FILTERED_GENOMES_FILE"]

include: "workflow/mash/Snakefile"
include: "workflow/anno_mlst_cdhit/Snakefile"
include: "workflow/infer_affinities/Snakefile"

# =============================================================================
# Default target: all notebook outputs
# =============================================================================

rule all:
    input:
        # Phase 1: Genome Filtering & Download
        f"{FIG}/1a_unfiltered_strains.png",
        f"{FIG}/1a_filtered_strains.png",
        f"{OUT}/data/1a_df_filtration.csv",
        f"{TEMP}/1b_genome_summary.csv",
        f"{TEMP}/1b_genome_metadata.csv",
        # Phase 2: Mash Filtration & Pangenome Construction
        f"{TEMP}/2a_genome_summary.csv",
        f"{OUT}/data/2b_mash_square.csv",
        f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
        f"{OUT}/data/2c_pangenome_summary.csv",
        f"{TEMP}/2d_enriched_metadata.csv",
        # Phase 3: Genome Characterization
        f"{OUT}/data/3a_car_summary.csv",
        f"{OUT}/data/3b_heaps_summary.csv",
        # Phase 4: NMF Decomposition
        f"{DATA}/processed/nmf-outputs/L.csv",
        f"{DATA}/processed/nmf-outputs/A.csv",
        f"{OUT}/data/4a_nmf_summary.csv",
        # Phase 5: Phylon Analysis
        f"{OUT}/data/5a_phylon_summary.csv",
        f"{OUT}/data/5b_gene_diff_stats.csv",
        f"{OUT}/data/5c_phylon_go_enrichment.csv",
        f"{OUT}/data/5d_genetic_variation.csv",
        f"{OUT}/data/5e_blast_summary.csv",
        # f"{OUT}/data/5f_affinity_summary.csv",  # excluded: no input strains in data/inferring_affinities/

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
# Phase 2: Mash Filtration & Pangenome Construction
# (Mash workflow rules included from workflow/mash/Snakefile)
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

checkpoint nb_2b_checkpoint:
    input:
        config="config.yml",
        summary=f"{TEMP}/2a_genome_summary.csv",
        metadata=f"{TEMP}/2a_genome_metadata.csv",
        mash_dist=f"{TEMP}/2b_mash/mash_distances.txt",
    output:
        summary=temp(f"{TEMP}/2b_genome_summary.csv"),
        metadata=temp(f"{TEMP}/2b_genome_metadata.csv"),
        mash_square=f"{OUT}/data/2b_mash_square.csv",
        mash_corr=f"{OUT}/data/2b_mash_corr_dist.csv",
        heatmap=f"{FIG}/2b_mash_heatmap.png",
        dist_all=f"{FIG}/2b_mash_dist_all.png",
        dist_filtered=f"{FIG}/2b_mash_dist_filtered.png",
        sensitivity=f"{FIG}/2b_sensitivity_analysis.png",
        clustermap_init=f"{FIG}/2b_clustermap_initial.png",
        clustermap_final=f"{FIG}/2b_clustermap_final.png",
        sizes_init=f"{FIG}/2b_cluster_sizes_initial.png",
        sizes_final=f"{FIG}/2b_cluster_sizes_final.png",
    shell:
        "python notebooks/2b_mash_filtration.py -- --config {input.config}"

rule nb_2c:
    input:
        config="config.yml",
        metadata=f"{TEMP}/2b_genome_metadata.csv",
        cdhit=f"{DATA}/processed/cd-hit-results/{SPECIES}",
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
# Phase 3: Genome Characterization
# (Anno/MLST/CD-HIT workflow rules included from workflow/anno_mlst_cdhit/Snakefile)
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

rule nb_4a:
    input:
        config="config.yml",
        acc=f"{DATA}/processed/CAR_genomes/df_acc.csv",
        metadata=f"{TEMP}/2d_enriched_metadata.csv",
    output:
        f"{DATA}/processed/nmf-outputs/L.csv",
        f"{DATA}/processed/nmf-outputs/A.csv",
        f"{FIG}/4a_mca_variance.png",
        f"{FIG}/4a_consensus_clustermap.png",
        f"{FIG}/4a_consensus_clustermap_filtered.png",
        f"{OUT}/data/4a_nmf_summary.csv",
        f"{OUT}/data/4a_consensus_clusters.csv",
    shell:
        "python notebooks/4a_nmf_decomposition.py -- --config {input.config}"

# =============================================================================
# Phase 5: Phylon Analysis
# (Infer affinities workflow rules included from workflow/infer_affinities/Snakefile)
# =============================================================================

rule nb_5a:
    input:
        config="config.yml",
        L=f"{DATA}/processed/nmf-outputs/L.csv",
        A=f"{DATA}/processed/nmf-outputs/A.csv",
        acc=f"{DATA}/processed/CAR_genomes/df_acc.csv",
        metadata=f"{TEMP}/2d_enriched_metadata.csv",
    output:
        f"{DATA}/processed/nmf-outputs/L_norm.csv",
        f"{DATA}/processed/nmf-outputs/A_norm.csv",
        f"{DATA}/processed/nmf-outputs/L_binarized.csv",
        f"{DATA}/processed/nmf-outputs/A_binarized.csv",
        f"{FIG}/5a_regression.png",
        f"{FIG}/5a_L_binarized_sorted.png",
        f"{FIG}/5a_A_binarized_sorted.png",
        f"{OUT}/data/5a_phylon_summary.csv",
    shell:
        "python notebooks/5a_phylon_characterization.py -- --config {input.config}"

rule nb_5b:
    input:
        config="config.yml",
        L_binarized=f"{DATA}/processed/nmf-outputs/L_binarized.csv",
    output:
        f"{FIG}/5b_clustermap.png",
        f"{FIG}/5b_phylon_dendrogram.png",
        f"{FIG}/5b_first_split_heatmap.png",
        f"{OUT}/data/5b_gene_diff_stats.csv",
    shell:
        "python notebooks/5b_gene_diff.py -- --config {input.config}"

rule nb_5c:
    input:
        config="config.yml",
        L_binarized=f"{DATA}/processed/nmf-outputs/L_binarized.csv",
    output:
        f"{DATA}/processed/all_functions.csv",
        f"{OUT}/data/5c_phylon_go_enrichment.csv",
        f"{FIG}/5c_enrichment_heatmap.png",
    shell:
        "python notebooks/5c_functional_enrichments.py -- --config {input.config}"

rule nb_5d:
    input:
        config="config.yml",
        L_binarized=f"{DATA}/processed/nmf-outputs/L_binarized.csv",
        A_binarized=f"{DATA}/processed/nmf-outputs/A_binarized.csv",
        pickle=f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
        metadata=f"{TEMP}/2d_enriched_metadata.csv",
    output:
        f"{FIG}/5d_gene_length_distribution.png",
        f"{OUT}/data/5d_genetic_variation.csv",
        f"{OUT}/data/5d_gene_alignment_summary.csv",
        f"{OUT}/data/5d_unique_genes_by_phylon.csv",
    shell:
        "python notebooks/5d_gene_alignment.py -- --config {input.config}"

rule nb_5e:
    input:
        config="config.yml",
        pickle=f"{DATA}/processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz",
    output:
        f"{OUT}/data/5e_blast_summary.csv",
    shell:
        "python notebooks/5e_blast_enrichment.py -- --config {input.config}"

rule nb_5f:
    input:
        config="config.yml",
        l_norm=f"{DATA}/processed/nmf-outputs/L_norm.csv",
        a_norm=f"{DATA}/processed/nmf-outputs/A_norm.csv",
        metadata=f"{TEMP}/2d_enriched_metadata.csv",
        p_new=f"{DATA}/inferring_affinities/combined_P_matrix.csv",
        mash_dist=f"{DATA}/inferring_affinities/combined_mash_distances.csv",
    output:
        f"{OUT}/data/5f_affinity_summary.csv",
    shell:
        "python notebooks/5f_infer_affinities.py -- --config {input.config}"
