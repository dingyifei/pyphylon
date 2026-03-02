import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import pandas as pd
    import yaml


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

    return CONFIG, DATA, TEMP


@app.cell
def _(DATA, TEMP):
    """Load 2b metadata and MLST report."""
    metadata_2b = pd.read_csv(
        os.path.join(TEMP, "2b_genome_metadata.csv"),
        dtype={"genome_id": str},
    )

    mlst_path = os.path.join(DATA, "processed", "mlst_report.txt")
    mlst_df = pd.read_csv(
        mlst_path,
        sep="\t",
        header=None,
        dtype="object",
        names=[
            "genome_id",
            "schema",
            "mlst",
            "allele1",
            "allele2",
            "allele3",
            "allele4",
            "allele5",
            "allele6",
            "allele7",
        ],
    )
    mlst_df["genome_id"] = mlst_df["genome_id"].apply(lambda x: os.path.basename(x).replace(".fna", ""))

    mo.output.replace(
        mo.md(
            f"Loaded inputs:\n\n"
            f"- **2b metadata:** {metadata_2b.shape[0]} genomes, {metadata_2b.shape[1]} columns\n"
            f"- **MLST report:** {mlst_df.shape[0]} entries"
        )
    )
    return metadata_2b, mlst_df


@app.cell
def _(metadata_2b, mlst_df):
    """Enrich metadata with MLST sequence types (vectorized merge)."""
    mlst_map = mlst_df.set_index("genome_id")["mlst"].replace("-", -1)
    enriched_metadata = metadata_2b.copy()
    enriched_metadata["mlst"] = enriched_metadata["genome_id"].map(mlst_map)

    _n_matched = enriched_metadata["mlst"].notna().sum()
    _mlst_counts = enriched_metadata["mlst"].value_counts().head(10)

    mo.output.replace(
        mo.md(
            f"MLST enrichment:\n\n"
            f"- Matched: **{_n_matched}** / {enriched_metadata.shape[0]} genomes\n"
            f"- Unique STs: **{enriched_metadata['mlst'].nunique()}**\n\n"
            f"Top 10 sequence types:\n\n" + _mlst_counts.to_markdown()
        )
    )
    return (enriched_metadata,)


@app.cell
def _(TEMP, enriched_metadata):
    """Save enriched metadata."""
    _output_path = os.path.join(TEMP, "2d_enriched_metadata.csv")
    enriched_metadata.to_csv(_output_path, index=False)

    mo.output.replace(
        mo.md(
            f"Saved `2d_enriched_metadata.csv`: **{enriched_metadata.shape[0]}** genomes, **{enriched_metadata.shape[1]}** columns"
        )
    )


if __name__ == "__main__":
    app.run()
