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
    mo.md(
        """
        # 2a: Clean Metadata

        De-duplicate genome metadata in preparation for **Mash** filtration
        and clustering. Duplicate entries sharing the same `biosample_accession`
        are removed so that each biological sample is represented exactly once.
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

    TEMP = CONFIG["TEMP_DIR"]

    os.makedirs(TEMP, exist_ok=True)

    return CONFIG, TEMP


@app.cell
def _(TEMP):
    """Load 1b outputs: genome summary and metadata."""
    input_summary = pd.read_csv(
        os.path.join(TEMP, "1b_genome_summary.csv"),
        dtype={"genome_id": str},
    )
    input_metadata = pd.read_csv(
        os.path.join(TEMP, "1b_genome_metadata.csv"),
        dtype={"genome_id": str},
    )

    mo.output.replace(
        mo.md(
            f"Loaded from 1b:\n\n"
            f"- **genome_summary:** {input_summary.shape[0]} rows\n"
            f"- **genome_metadata:** {input_metadata.shape[0]} rows"
        )
    )
    return input_metadata, input_summary


@app.cell
def _():
    mo.md(
        """
        ## De-duplicate Entries

        Ensure `biosample_accession` is unique — drop rows that share the same
        biosample so downstream analyses are not biased by redundant sequences.
        The genome summary table is then synced to match.
        """
    )


@app.cell
def _(input_metadata, input_summary):
    """De-duplicate metadata by biosample_accession and sync summary."""
    dedup_metadata = input_metadata.drop_duplicates(subset=["biosample_accession"])
    dedup_summary = input_summary[input_summary["genome_id"].isin(dedup_metadata["genome_id"])]

    n_removed = len(input_metadata) - len(dedup_metadata)

    mo.output.replace(
        mo.md(
            f"De-duplication by `biosample_accession`:\n\n"
            f"- Removed **{n_removed}** duplicate rows\n"
            f"- **dedup_summary:** {dedup_summary.shape[0]} rows\n"
            f"- **dedup_metadata:** {dedup_metadata.shape[0]} rows"
        )
    )
    return dedup_metadata, dedup_summary


@app.cell
def _():
    mo.md("## Save Files")


@app.cell
def _(TEMP, dedup_metadata, dedup_summary):
    """Save de-duplicated CSVs for downstream notebooks."""
    dedup_summary.to_csv(os.path.join(TEMP, "2a_genome_summary.csv"), index=False)
    dedup_metadata.to_csv(os.path.join(TEMP, "2a_genome_metadata.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"Saved:\n\n"
            f"- `{TEMP}2a_genome_summary.csv` ({dedup_summary.shape[0]} rows)\n"
            f"- `{TEMP}2a_genome_metadata.csv` ({dedup_metadata.shape[0]} rows)"
        )
    )
    return


if __name__ == "__main__":
    app.run()
