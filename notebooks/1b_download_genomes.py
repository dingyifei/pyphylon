import marimo

__generated_with = "0.20.2"
app = marimo.App(width="full")

with app.setup:
    import os
    import sys

    import marimo as mo
    import pandas as pd
    import yaml
    from tqdm.auto import tqdm

    from pyphylon.downloads import download_genomes_bvbrc
    from pyphylon.util import remove_empty_files


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
    DEBUG = CONFIG["DEBUG"]

    os.makedirs(TEMP, exist_ok=True)
    os.makedirs(os.path.join(DATA, "raw", "genomes"), exist_ok=True)
    return DATA, TEMP


@app.cell
def _(TEMP):
    """Load 1a outputs: genome summary and metadata."""
    input_summary = pd.read_csv(
        os.path.join(TEMP, "1a_genome_summary.csv"),
        dtype={"genome_id": str},
    )
    input_metadata = pd.read_csv(
        os.path.join(TEMP, "1a_genome_metadata.csv"),
        dtype={"genome_id": str},
    )

    mo.output.replace(
        mo.md(
            f"Loaded from 1a:\n\n"
            f"- **genome_summary:** {input_summary.shape[0]} rows\n"
            f"- **genome_metadata:** {input_metadata.shape[0]} rows"
        )
    )
    return input_metadata, input_summary


@app.cell
def _(DATA, input_summary):
    """Download genome .fna files from BV-BRC. Skips already-downloaded files."""
    RAW_GENOMES = os.path.join(DATA, "raw", "genomes")

    bad_genomes = download_genomes_bvbrc(
        genomes=input_summary["genome_id"],
        output_dir=RAW_GENOMES,
        filetypes=["fna"],
    )

    mo.output.replace(mo.md(f"Download complete. **{len(bad_genomes)}** genomes failed to download."))
    return RAW_GENOMES, bad_genomes


@app.cell
def _(RAW_GENOMES):
    """Remove empty files from download subdirectories."""
    empty_files = []
    for subdir in tqdm(os.listdir(RAW_GENOMES), desc="Checking for empty files"):
        subdir_path = os.path.join(RAW_GENOMES, subdir)
        if os.path.isdir(subdir_path):
            files = remove_empty_files(subdir_path)
            empty_files.extend(files)

    mo.output.replace(mo.md(f"Removed **{len(empty_files)}** empty files."))
    return


@app.cell
def _(bad_genomes, input_metadata, input_summary):
    """Filter summary and metadata to keep only successfully downloaded genomes."""
    downloaded_genomes = set(input_summary.genome_id.astype(str)) - set(bad_genomes)

    dl_summary = (
        input_summary.drop_duplicates(subset=["genome_id"])
        .set_index("genome_id")
        .loc[sorted(downloaded_genomes)]
        .reset_index()
    )

    dl_metadata = (
        input_metadata.drop_duplicates(subset=["genome_id"])
        .set_index("genome_id")
        .loc[sorted(downloaded_genomes)]
        .reset_index()
    )

    mo.output.replace(
        mo.md(
            f"After filtering to downloaded genomes:\n\n"
            f"- **dl_summary:** {dl_summary.shape[0]} rows\n"
            f"- **dl_metadata:** {dl_metadata.shape[0]} rows"
        )
    )
    return dl_metadata, dl_summary


@app.cell
def _(TEMP, dl_metadata, dl_summary):
    """Save filtered CSVs for downstream notebooks."""
    dl_summary.to_csv(os.path.join(TEMP, "1b_genome_summary.csv"), index=False)
    dl_metadata.to_csv(os.path.join(TEMP, "1b_genome_metadata.csv"), index=False)

    mo.output.replace(
        mo.md(
            f"Saved:\n\n"
            f"- `{TEMP}1b_genome_summary.csv` ({dl_summary.shape[0]} rows)\n"
            f"- `{TEMP}1b_genome_metadata.csv` ({dl_metadata.shape[0]} rows)"
        )
    )
    return


if __name__ == "__main__":
    app.run()
