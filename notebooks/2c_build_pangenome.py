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

    from pyphylon.pangenome import build_cds_pangenome

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

    TEMP = CONFIG["TEMP_DIR"]
    DATA = CONFIG["DATA_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]
    SPECIES = CONFIG["PG_NAME"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)

    return CONFIG, DATA, FIG, OUT, SPECIES, TEMP


@app.cell
def _(TEMP):
    """Load mash-filtered genome metadata from 2b."""
    metadata_2b = pd.read_csv(
        os.path.join(TEMP, "2b_genome_metadata.csv"),
        dtype={"genome_id": str},
    )

    mo.output.replace(
        mo.md(
            f"Loaded `2b_genome_metadata.csv`: **{metadata_2b.shape[0]}** genomes, **{metadata_2b.shape[1]}** columns"
        )
    )
    return (metadata_2b,)


@app.cell
def _(DATA, metadata_2b):
    """Discover BAKTA .faa files and filter to mash-scrubbed genomes."""
    bakta_dir = os.path.join(DATA, "processed", "bakta")

    # Enumerate all .faa files from BAKTA output directories
    bakta_faa_paths = []
    if os.path.isdir(bakta_dir):
        bakta_faa_paths = [
            os.path.join(bakta_dir, folder, folder + ".faa")
            for folder in os.listdir(bakta_dir)
            if os.path.isdir(os.path.join(bakta_dir, folder))
        ]

    # Filter to genomes that passed mash filtration
    genome_id_set = set(metadata_2b["genome_id"].astype(str))
    filtered_faa_paths = [
        f for f in bakta_faa_paths if os.path.basename(os.path.dirname(f)) in genome_id_set and os.path.isfile(f)
    ]

    mo.output.replace(
        mo.md(
            f"BAKTA .faa discovery:\n\n"
            f"- Total BAKTA directories: **{len(bakta_faa_paths)}**\n"
            f"- Matched to 2b metadata: **{len(filtered_faa_paths)}** / {metadata_2b.shape[0]} genomes"
        )
    )
    return (filtered_faa_paths,)


@app.cell
def _(CONFIG, DATA, SPECIES, filtered_faa_paths):
    """Build CDS pangenome using CD-HIT (or load existing results)."""
    cdhit_output_dir = os.path.join(DATA, "processed", "cd-hit-results")
    os.makedirs(cdhit_output_dir, exist_ok=True)

    pickle_genes = os.path.join(cdhit_output_dir, f"{SPECIES}_strain_by_gene.pickle.gz")
    pickle_alleles = os.path.join(cdhit_output_dir, f"{SPECIES}_strain_by_allele.pickle.gz")

    if os.path.exists(pickle_genes) and CONFIG.get("REUSE_TEMP", True):
        # Load pre-existing results
        df_genes = pd.read_pickle(pickle_genes)
        df_alleles = pd.read_pickle(pickle_alleles)
        mo.output.replace(
            mo.md(
                f"Loaded existing pangenome from pickle:\n\n"
                f"- **Gene matrix:** {df_genes.shape[0]} genes x {df_genes.shape[1]} genomes\n"
                f"- **Allele matrix:** {df_alleles.shape[0]} alleles x {df_alleles.shape[1]} genomes"
            )
        )
    elif len(filtered_faa_paths) == 0:
        mo.output.replace(
            mo.md(
                "**No BAKTA .faa files found.** Run the anno_mlst_cdhit workflow first.\n\n"
                "```bash\n"
                "cd workflow/anno_mlst_cdhit\n"
                "snakemake -d ../../data --use-singularity -c 10\n"
                "```"
            )
        )
        df_genes = pd.DataFrame()
        df_alleles = pd.DataFrame()
    else:
        # Run CD-HIT pangenome construction
        df_alleles, df_genes, _header_to_allele = build_cds_pangenome(
            genome_faa_paths=filtered_faa_paths,
            output_dir=cdhit_output_dir,
            name=SPECIES,
            cdhit_args={"-n": 5, "-c": 0.8, "-aL": 0.8, "-T": 0, "-M": 0},
            fastasort_path=None,
            save_csv=False,
        )
        mo.output.replace(
            mo.md(
                f"Built CDS pangenome with CD-HIT:\n\n"
                f"- **Gene matrix:** {df_genes.shape[0]} genes x {df_genes.shape[1]} genomes\n"
                f"- **Allele matrix:** {df_alleles.shape[0]} alleles x {df_alleles.shape[1]} genomes"
            )
        )

    return df_alleles, df_genes


@app.cell
def _(FIG, df_genes):
    """Plot genes-per-genome distribution."""
    if df_genes.empty:
        mo.output.replace(mo.md("*No pangenome data to plot.*"))
    else:
        _genes_per_genome = df_genes.fillna(0).sum(axis=0).astype(int)

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(_genes_per_genome, bins=40, edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Number of genes")
        ax.set_ylabel("Number of genomes")
        ax.set_title("Genes per Genome Distribution")
        ax.axvline(_genes_per_genome.mean(), color="red", linestyle="--", label=f"Mean: {_genes_per_genome.mean():.0f}")
        ax.legend()
        plt.tight_layout()

        fig.savefig(os.path.join(FIG, "2c_genes_per_genome.png"), bbox_inches="tight")
        mo.output.replace(fig)


@app.cell
def _(OUT, df_alleles, df_genes):
    """Save pangenome summary statistics."""
    if df_genes.empty:
        mo.output.replace(mo.md("*No pangenome data to summarize.*"))
    else:
        _genes_per_genome = df_genes.fillna(0).sum(axis=0).astype(int)
        _summary = pd.DataFrame(
            {
                "metric": [
                    "n_genomes",
                    "n_genes",
                    "n_alleles",
                    "mean_genes_per_genome",
                    "std_genes_per_genome",
                    "min_genes_per_genome",
                    "max_genes_per_genome",
                ],
                "value": [
                    df_genes.shape[1],
                    df_genes.shape[0],
                    df_alleles.shape[0],
                    f"{_genes_per_genome.mean():.1f}",
                    f"{_genes_per_genome.std():.1f}",
                    f"{_genes_per_genome.min():.0f}",
                    f"{_genes_per_genome.max():.0f}",
                ],
            }
        )

        _summary.to_csv(os.path.join(OUT, "data", "2c_pangenome_summary.csv"), index=False)

        mo.output.replace(
            mo.md(
                "Pangenome summary saved to `output/data/2c_pangenome_summary.csv`:\n\n"
                "| Metric | Value |\n"
                "|--------|-------|\n" + "\n".join(f"| {row.metric} | {row.value} |" for row in _summary.itertuples())
            )
        )


if __name__ == "__main__":
    app.run()
