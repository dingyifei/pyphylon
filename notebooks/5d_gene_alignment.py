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
    from tqdm.auto import tqdm

    from pyphylon.biointerp import get_pg_to_locus_map
    from pyphylon.plotting import plot_circular_genome
    from pyphylon.plotting_util import (
        adjust_gene_order,
        check_strict_sequence,
        create_gene_count_between_anchor_genes_for_all,
        create_strain_groups,
        find_once_genes,
        generate_gene_names,
        get_reference_order,
        gff2pandas,
        identify_genetic_variation,
        reorder_to_start_with_one,
        unique_genes_by_phylon,
        update_strain_vector,
    )

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    mo.md(
        """
        # 5d: Gene Alignment Visualization

        Align genes across Complete genomes using **anchor genes** — genes appearing
        exactly once in every strain. Strains are grouped by gene order consistency,
        reversed if needed, and rotated to start at gene 1. This standardized ordering
        enables circular genome visualization of phylon-specific gene locations.

        **Pipeline:** Parse GFF annotations → identify anchor genes → standardize gene
        order → classify genetic variation → generate circular genome plots per phylon.
        """
    )


@app.cell
def _():
    mo.md("## Setup")


@app.cell
def _():
    """Parse config and set up directories."""
    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    DATA = CONFIG["DATA_DIR"]
    TEMP = CONFIG["TEMP_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]
    SPECIES = CONFIG["PG_NAME"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)

    return CONFIG, DATA, FIG, OUT, SPECIES, TEMP


@app.cell
def _():
    mo.md(
        r"""
        ## Load Inputs

        Load the P matrix, enriched metadata (filtered to Complete genomes), L/A
        binarized matrices, and the pangenome cluster → locus\_tag mapping.
        """
    )


@app.cell
def _(DATA, SPECIES, TEMP):
    """Load P matrix (for genome filtering), metadata, L/A binarized, and pg2locus_map."""
    # P matrix — only used to filter metadata to genomes present in the pangenome
    df_genes = pd.read_pickle(os.path.join(DATA, f"processed/cd-hit-results/{SPECIES}_strain_by_gene.pickle.gz"))

    # Enriched metadata, filtered to Complete genomes present in P matrix
    metadata_raw = pd.read_csv(os.path.join(TEMP, "2d_enriched_metadata.csv"), dtype={"genome_id": str})
    metadata_complete = metadata_raw[
        metadata_raw["genome_id"].isin(df_genes.columns) & (metadata_raw["genome_status"] == "Complete")
    ].copy()

    # Auto-select reference strain (first Complete genome)
    REFERENCE_STRAIN = metadata_complete["genome_id"].iloc[0]

    # L and A binarized matrices
    L_binarized_clean = pd.read_csv(os.path.join(DATA, "processed/nmf-outputs/L_binarized.csv"), index_col=0)
    # Defensive index cleanup: strip 'yogenes' prefix (no-op for CJejuni, needed for SPyogenes)
    L_binarized_clean.index = L_binarized_clean.index.str.replace("yogenes", "", regex=False)

    A_binarized = pd.read_csv(os.path.join(DATA, "processed/nmf-outputs/A_binarized.csv"), index_col=0)

    # Pangenome cluster → locus_tag mapping (from allele_names.tsv)
    pg2locus_map = get_pg_to_locus_map(DATA, SPECIES)

    n_complete = len(metadata_complete)
    mo.output.replace(
        mo.md(
            f"Loaded inputs:\n\n"
            f"- **Complete genomes in P matrix:** {n_complete}\n"
            f"- **Reference strain:** `{REFERENCE_STRAIN}`\n"
            f"- **L_binarized:** {L_binarized_clean.shape[0]} genes x {L_binarized_clean.shape[1]} phylons\n"
            f"- **A_binarized:** {A_binarized.shape[0]} phylons x {A_binarized.shape[1]} strains\n"
            f"- **pg2locus_map:** {len(pg2locus_map)} mappings"
        )
    )
    return A_binarized, L_binarized_clean, REFERENCE_STRAIN, metadata_complete, pg2locus_map


@app.cell
def _():
    mo.md(
        """
        ## Parse GFF Files

        Parse BAKTA GFF3 annotations for each Complete genome. For each strain, merge
        gene positions with pangenome cluster IDs and keep only the main chromosome
        (largest accession). The result is an ordered gene vector per strain.
        """
    )


@app.cell
def _(DATA, metadata_complete, pg2locus_map):
    """Parse GFF files for all Complete genomes and build strain_vectors (cached)."""
    with mo.persistent_cache("5d_gff_parsing"):
        strain_vectors = {}
        for strain in tqdm(metadata_complete["genome_id"], desc="Parsing GFF files"):
            gff_path = os.path.join(DATA, f"processed/bakta/{strain}/{strain}.gff3")
            DF_gff, _size, _oric = gff2pandas(gff_path)
            DF_gff = pd.merge(DF_gff, pg2locus_map, left_on="locus_tag", right_on="gene_id", how="left")
            DF_gff.rename(columns={"cluster": "gene"}, inplace=True)
            # Keep only the largest accession (main chromosome)
            DF_gff = DF_gff[DF_gff["accession"] == DF_gff["accession"].value_counts().index[0]]
            DF_gff = DF_gff[["gene", "start"]]
            strain_vectors[strain] = DF_gff.sort_values("start")["gene"].tolist()

    mo.output.replace(
        mo.md(
            f"Parsed GFF files for **{len(strain_vectors)}** strains\n\n"
            f"- Median genes per strain: {sorted(len(g) for g in strain_vectors.values())[len(strain_vectors) // 2]}"
        )
    )
    return (strain_vectors,)


@app.cell
def _():
    mo.md(
        """
        ## Gene Order Standardization

        Identify anchor genes (single-copy genes shared by all strains), group strains
        by consistent gene ordering, establish a canonical reference order from the
        largest group, then renumber, reverse, and rotate each strain's gene vector
        so they all start at gene 1.
        """
    )


@app.cell
def _(REFERENCE_STRAIN, strain_vectors):
    """Standardize gene order: filter, find once-genes, group strains, build reference order, reorder."""
    # Filter to strains with <= 6000 genes
    strain_vectors_filtered = {strain: genes for strain, genes in strain_vectors.items() if len(genes) <= 6000}

    # Find genes appearing exactly once in every strain (anchor genes)
    common_gene_count, once_gene_count, once_genes = find_once_genes(strain_vectors_filtered)

    # Group strains by consistent gene ordering relative to reference
    strain_groups = create_strain_groups(strain_vectors_filtered, once_genes, REFERENCE_STRAIN)

    # Use the largest group as reference for gene ordering
    largest_group = max(strain_groups, key=lambda k: len(strain_groups[k]))
    strain_vectors_reference = {
        k: strain_vectors_filtered[k] for k in strain_groups[largest_group] if k in strain_vectors_filtered
    }

    # Get canonical gene order from reference group
    reference_ordered_genes = get_reference_order(strain_vectors_reference, once_genes)

    # Apply reference mapping: replace gene names with numeric positions
    updated_strain_vectors = update_strain_vector(reference_ordered_genes, strain_vectors_filtered)

    # Reverse generally-decreasing gene lists
    strain_vectors_reordered, count_reversed = adjust_gene_order(updated_strain_vectors)

    # Rotate gene lists to start with gene '1'
    strain_vectors_final, count_reordered = reorder_to_start_with_one(strain_vectors_reordered)

    # Validate strict sequence
    _seq_results, total_true, total_false = check_strict_sequence(strain_vectors_final)

    mo.output.replace(
        mo.md(
            f"Gene order standardization:\n\n"
            f"- **Filtered strains (<=6000 genes):** {len(strain_vectors_filtered)}\n"
            f"- **Common genes:** {common_gene_count}\n"
            f"- **Once-genes (anchor genes):** {once_gene_count}\n"
            f"- **Strain groups:** {len(strain_groups)} "
            f"(largest: {largest_group} with {len(strain_groups[largest_group])} strains)\n"
            f"- **Reference ordered genes:** {len(reference_ordered_genes)}\n"
            f"- **Strains reversed:** {count_reversed}\n"
            f"- **Strains re-rotated:** {count_reordered}\n"
            f"- **Strictly ordered:** {total_true} | Other orders: {total_false}"
        )
    )
    return once_genes, strain_vectors_final


@app.cell
def _():
    mo.md(
        """
        ## Genetic Variation

        Assign descriptive names to non-anchor genes based on their flanking anchor
        positions, count genes between each anchor pair, and classify each strain's
        structural variation type: no variation, inversion, translocation, or other.
        """
    )


@app.cell
def _(strain_vectors_final):
    """Generate gene names relative to anchors and identify genetic variation types."""
    # Descriptive gene names relative to anchor genes (e.g. "42_3_2_43")
    gene_mapping_to_anchor_genes = generate_gene_names(strain_vectors_final)

    # Count genes between each pair of anchor genes for all strains
    gene_count_between_anchors = create_gene_count_between_anchor_genes_for_all(gene_mapping_to_anchor_genes)

    # Identify variation types (no variation, inversion, translocation, others) per strain
    variation_df = identify_genetic_variation(strain_vectors_final)

    var_counts = variation_df["Variation"].value_counts()
    mo.output.replace(
        mo.md(
            f"Gene naming and variation analysis:\n\n"
            f"- **Gene mapping matrix:** {gene_mapping_to_anchor_genes.shape}\n"
            f"- **Strains with anchor gene counts:** {len(gene_count_between_anchors)}\n"
            f"- **Variation types:**\n" + "\n".join(f"  - {vtype}: {count}" for vtype, count in var_counts.items())
        )
    )
    return gene_mapping_to_anchor_genes, variation_df


@app.cell
def _():
    mo.md(
        """
        ## Gene Length Distribution

        Histogram of gene counts per Complete genome after GFF parsing.
        """
    )


@app.cell
def _(FIG, strain_vectors):
    """Plot histogram of gene counts per strain."""
    gene_lengths = [len(genes) for genes in strain_vectors.values()]

    fig_hist, ax_hist = plt.subplots(figsize=(8, 5))
    ax_hist.hist(gene_lengths, bins=30, color="steelblue", edgecolor="black")
    ax_hist.set_title("Distribution of Gene Counts per Complete Genome")
    ax_hist.set_xlabel("Gene Count")
    ax_hist.set_ylabel("Frequency")
    fig_hist.tight_layout()

    fig_hist.savefig(os.path.join(FIG, "5d_gene_length_distribution.png"), bbox_inches="tight")
    plt.close(fig_hist)

    mo.output.replace(mo.md(f"Saved `5d_gene_length_distribution.png` ({len(gene_lengths)} strains)"))


@app.cell
def _():
    mo.md(
        """
        ## Circular Genome Plots

        Circular plots show gene positions for each phylon. Orange/coral segments
        indicate phylon-specific genes; blue segments indicate anchor genes (conserved
        single-copy genes used for ordering).
        """
    )


@app.cell
def _(A_binarized, FIG, L_binarized_clean, strain_vectors_final):
    """Generate circular genome plots for each phylon and compute unique genes per phylon."""
    unique_genes_dict = unique_genes_by_phylon(L_binarized_clean)

    phylon_cols = L_binarized_clean.columns.tolist()
    phylon_tabs = {}

    for phylon_col in phylon_cols:
        # Get genes belonging to this phylon
        phylon_genes = L_binarized_clean.index[L_binarized_clean[phylon_col] == 1].tolist()
        if len(phylon_genes) == 0:
            continue

        # Find a strain belonging to this phylon (from A_binarized)
        plot_strain = None
        if phylon_col in A_binarized.index:
            phylon_strains = A_binarized.columns[A_binarized.loc[phylon_col] == 1].tolist()
            for s in phylon_strains:
                if s in strain_vectors_final:
                    plot_strain = s
                    break

        # Fall back to first available strain
        if plot_strain is None:
            plot_strain = next(iter(strain_vectors_final))

        # Pass a COPY of the gene list (plot_circular_genome rotates it in-place)
        fig_circ = plot_circular_genome(list(strain_vectors_final[plot_strain]), phylon_genes, phylon_col, plot_strain)
        fig_circ.write_image(os.path.join(FIG, f"5d_circular_{phylon_col}.png"))
        phylon_tabs[phylon_col] = mo.ui.plotly(fig_circ)

    mo.output.replace(
        mo.vstack(
            [
                mo.md(f"Saved **{len(phylon_tabs)}** circular genome plots to `output/figures/5d_circular_*.png`"),
                mo.ui.tabs(phylon_tabs),
            ]
        )
    )
    return (unique_genes_dict,)


@app.cell
def _():
    mo.md(
        """
        ## Unique Genes by Phylon

        Genes that belong exclusively to a single phylon (not shared with any other).
        """
    )


@app.cell
def _(unique_genes_dict):
    """Display unique gene counts per phylon as an interactive table."""
    _rows = [
        {"phylon": p, "unique_genes": len(genes)}
        for p, genes in sorted(unique_genes_dict.items(), key=lambda x: len(x[1]), reverse=True)
    ]
    mo.output.replace(mo.ui.table(pd.DataFrame(_rows), label="Unique gene counts by phylon"))


@app.cell
def _():
    mo.md("## Save Results")


@app.cell
def _(OUT, strain_vectors_final, unique_genes_dict, variation_df):
    """Save all data outputs and display final summary."""
    # Genetic variation per strain
    variation_df.to_csv(os.path.join(OUT, "data", "5d_genetic_variation.csv"), index=False)

    # Summary metrics
    n_anchor = 0
    for genes in strain_vectors_final.values():
        n_anchor = sum(1 for g in genes if isinstance(g, int))
        break  # All strains have the same anchor count

    var_counts = variation_df["Variation"].value_counts()
    summary = pd.DataFrame(
        {
            "metric": [
                "n_strains_final",
                "n_anchor_genes",
                "n_no_variation",
                "n_inversion",
                "n_translocation",
                "n_other",
            ],
            "value": [
                len(strain_vectors_final),
                n_anchor,
                int(var_counts.get("no variation", 0)),
                int(var_counts.get("inversion", 0)),
                int(var_counts.get("translocation", 0)),
                int(var_counts.get("others", 0)),
            ],
        }
    )
    summary.to_csv(os.path.join(OUT, "data", "5d_gene_alignment_summary.csv"), index=False)

    # Unique genes per phylon
    unique_genes_rows = []
    for phylon, genes in unique_genes_dict.items():
        for gene in genes:
            unique_genes_rows.append({"phylon": phylon, "gene": gene})
    if unique_genes_rows:
        pd.DataFrame(unique_genes_rows).to_csv(os.path.join(OUT, "data", "5d_unique_genes_by_phylon.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"**5d summary:**\n\n"
            f"- Strains with standardized gene order: {len(strain_vectors_final)}\n"
            f"- Anchor genes: {n_anchor}\n"
            f"- Variation: "
            + ", ".join(f"{v}: {c}" for v, c in var_counts.items())
            + f"\n- Unique gene-phylon pairs: {len(unique_genes_rows)}\n\n"
            f"**Outputs:**\n"
            f"- `5d_gene_length_distribution.png`\n"
            f"- `5d_circular_*.png` (one per phylon)\n"
            f"- `5d_genetic_variation.csv`\n"
            f"- `5d_gene_alignment_summary.csv`\n"
            f"- `5d_unique_genes_by_phylon.csv`"
        )
    )


if __name__ == "__main__":
    app.run()
