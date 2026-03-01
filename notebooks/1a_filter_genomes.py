import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import matplotlib.pyplot as plt
    import seaborn as sns
    import yaml

    import pyphylon.qcqa as qcqa
    from pyphylon.downloads import get_scaffold_n50_for_species, query_bvbrc_genomes


@app.cell
def _():
    """Parse config and set up directories."""
    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    TEMP = CONFIG["TEMP_DIR"]
    OUTPUT = CONFIG["OUTPUT_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    SPECIES_NAME = CONFIG["SPECIES_NAME"]
    TAXON_ID = CONFIG["TAXON_ID"]
    DEBUG = CONFIG["DEBUG"]

    os.makedirs(TEMP, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT, "data"), exist_ok=True)

    # Plot defaults
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")

    return CONFIG, DEBUG, FIG, OUTPUT, SPECIES_NAME, TAXON_ID, TEMP


@app.cell
def _(SPECIES_NAME, TAXON_ID):
    """Query BV-BRC API and filter by species. Cached because API calls are slow."""
    with mo.persistent_cache("1a_query_genomes"):
        genome_df = query_bvbrc_genomes(TAXON_ID, genome_status="Complete", genome_quality="Good")
        genome_df = genome_df.set_index("genome_id")
        species_summary = qcqa.filter_by_species(genome_df.copy(), SPECIES_NAME)
        scaffold_n50 = get_scaffold_n50_for_species(species_summary.taxon_id.mode().values[0])

    mo.output.replace(
        mo.md(
            f"Retrieved **{genome_df.shape[0]}** genomes from BV-BRC API\n\n"
            f"Species-filtered: **{species_summary.shape[0]}** genomes\n\n"
            f"Reference scaffold N50: **{scaffold_n50:,}** bp"
        )
    )
    return genome_df, scaffold_n50, species_summary


@app.cell
def _(FIG, scaffold_n50, species_summary):
    """Plot unfiltered strain distributions and save figures."""
    # --- Unfiltered jointplot ---
    h_unfiltered = sns.jointplot(
        data=species_summary,
        x="genome_length",
        y="patric_cds",
        hue="genome_status",
        alpha=0.75,
        height=4,
    )
    h_unfiltered.ax_joint.legend(title="BV-BRC\nstrain type")
    h_unfiltered.ax_joint.set_xlabel("genome length")
    h_unfiltered.ax_joint.set_ylabel("BV-BRC predicted gene count")
    h_unfiltered.savefig(os.path.join(FIG, "1a_unfiltered_strains.png"), bbox_inches="tight")

    # --- N50 histogram for Complete genomes ---
    species_complete = species_summary[species_summary.genome_status == "Complete"]
    min_thresh_n50 = int(0.85 * scaffold_n50)

    fig_n50, ax_n50 = plt.subplots()
    sns.histplot(species_complete.contig_n50.dropna().astype("int"), ax=ax_n50)
    ax_n50.axvline(x=min_thresh_n50, color="#ff00ff", linestyle="--")
    fig_n50.savefig(os.path.join(FIG, "1a_unfiltered_n50.png"), bbox_inches="tight")

    # --- CheckM completeness & contamination line plots ---
    fig_checkm, axes = plt.subplots(1, 2, figsize=(10, 3.5))

    completeness = species_summary["checkm_completeness"].dropna().sort_values(ascending=False).reset_index(drop=True)
    axes[0].plot(completeness.values, linewidth=0.8)
    axes[0].set_xlabel("Genome (sorted)")
    axes[0].set_ylabel("CheckM Completeness (%)")
    axes[0].set_title("Completeness")
    axes[0].axhline(y=92, color="#ff00ff", linestyle="--", label="92% cutoff")
    axes[0].legend()

    contamination = species_summary["checkm_contamination"].dropna().sort_values().reset_index(drop=True)
    axes[1].plot(contamination.values, linewidth=0.8, color="tab:orange")
    axes[1].set_xlabel("Genome (sorted)")
    axes[1].set_ylabel("CheckM Contamination (%)")
    axes[1].set_title("Contamination")
    axes[1].axhline(y=5, color="#ff00ff", linestyle="--", label="5% cutoff")
    axes[1].legend()
    fig_checkm.tight_layout()

    mo.output.replace(mo.vstack([h_unfiltered.figure, fig_n50, fig_checkm]))
    return (min_thresh_n50,)


@app.cell
def _(DEBUG, genome_df, min_thresh_n50, species_summary):
    """Apply genome quality filters. DEBUG truncation is in this cell to avoid redefinition."""
    filtered_species_summary, df_filtration = qcqa.filter_by_genome_quality(
        species_summary,
        min_thresh_n50=min_thresh_n50,
        max_contig=None,
        contamination_cutoff=5.0,
        completeness_cutoff=92.0,
        checkm_filter_statuses=("Complete", "WGS"),
        checkm_missing="drop",
        return_stats=True,
    )

    # DEBUG mode: truncate to 50 genomes for fast iteration
    if DEBUG:
        filtered_species_summary = filtered_species_summary[:50]

    # Derive metadata from the original unfiltered genome_df for filtered indices
    filtered_species_metadata = genome_df.loc[filtered_species_summary.index]

    mo.output.replace(
        mo.vstack(
            [
                mo.md(f"**Filtered strains:** {filtered_species_summary.shape[0]}"),
                mo.md("### Filtration Report"),
                mo.ui.table(df_filtration),
            ]
        )
    )
    return df_filtration, filtered_species_metadata, filtered_species_summary


@app.cell
def _(FIG, filtered_species_summary):
    """Plot filtered strain distributions and save figures."""
    # --- Filtered jointplot ---
    h_filtered = sns.jointplot(
        data=filtered_species_summary,
        x="genome_length",
        y="patric_cds",
        hue="genome_status",
        alpha=0.75,
        height=4,
    )
    h_filtered.ax_joint.legend(title="BV-BRC\nstrain type")
    h_filtered.ax_joint.set_xlabel("genome length")
    h_filtered.ax_joint.set_ylabel("BV-BRC predicted gene count")
    h_filtered.savefig(os.path.join(FIG, "1a_filtered_strains.png"), bbox_inches="tight")

    # --- GC content vs contigs (diagnostic, not saved) ---
    h_gc = sns.jointplot(
        data=filtered_species_summary,
        x="gc_content",
        y="contigs",
        hue="genome_status",
        alpha=0.75,
        height=4,
    )
    h_gc.ax_joint.legend(title="BV-BRC\nstrain type", bbox_to_anchor=(1.45, 1.4))
    h_gc.ax_joint.set_xlabel("GC Content")
    h_gc.ax_joint.set_ylabel("number of contigs")

    mo.output.replace(mo.vstack([h_filtered.figure, h_gc.figure]))
    return


@app.cell
def _(OUTPUT, TEMP, df_filtration, filtered_species_metadata, filtered_species_summary):
    """Save all CSV outputs."""
    filtered_species_summary.to_csv(os.path.join(TEMP, "1a_genome_summary.csv"))
    filtered_species_metadata.to_csv(os.path.join(TEMP, "1a_genome_metadata.csv"))
    df_filtration.to_csv(os.path.join(OUTPUT, "data", "1a_df_filtration.csv"))

    mo.output.replace(
        mo.md(
            f"Saved:\n"
            f"- `{TEMP}1a_genome_summary.csv` ({filtered_species_summary.shape[0]} rows)\n"
            f"- `{TEMP}1a_genome_metadata.csv` ({filtered_species_metadata.shape[0]} rows)\n"
            f"- `{OUTPUT}data/1a_df_filtration.csv`"
        )
    )
    return


if __name__ == "__main__":
    app.run()
