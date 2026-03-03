# Examples (Legacy)

This directory contains legacy configuration files and example data for the pyphylon pipeline. The original Jupyter notebooks have been converted to **marimo** notebooks in `notebooks/`.

See the [root README](../README.md) for the current pipeline workflow.

## Example Data (S. pyogenes)

Two options for obtaining example data:

1. **Metadata only** (run bioinformatics yourself):
   Download [metadata](https://www.dropbox.com/s/ksvdgi8xgfx5m2r/spyogenes_metadata_summary.tar?dl=0) and extract to `data/metadata/`

2. **Pre-computed results** (skip bioinformatics):
   Download [all data](https://www.dropbox.com/scl/fi/fam8tnk8hhfgek3fwnn9l/SPyogenes_example_v2.tar?rlkey=lom35fwuy8sd5y2181if3aiwo&st=rsu9z4fg&dl=0) and extract to `data/`

## Using Example Data with the New Pipeline

1. Copy extracted data into the root `data/` directory
2. Update `config.yml` with the appropriate species settings (see `config.example.yml`)
3. Run the pipeline from the repository root:
   ```bash
   # Host mode
   conda run -n pyphylon-marimo python notebooks/1a_filter_genomes.py -- --config config.yml

   # Docker mode (full pipeline)
   docker compose run pipeline snakemake --cores 4 --sdm apptainer
   ```
