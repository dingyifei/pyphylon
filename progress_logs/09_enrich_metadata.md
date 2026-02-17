# 09 — Update Notebook 2d for Current Project Structure

**Date**: 2026-02-17 10:09 AM PST

## Summary

Updated `examples/2d_enrich_metadata.ipynb` to use the `temp/` and `output/` directory convention, matching notebooks 1a through 2c.

## Key Steps

1. **Cell 1 (config setup)** — Added `temp_folder` and `data_dir` variables derived from `config.yml`, matching the pattern from notebook 2c.
2. **Cell 2 (path definitions)** — Fixed `MLST` path to use `data_dir` (`data/processed/mlst_report.txt`) and `METADATA` path to use `temp_folder` (`../temp/2b_genome_metadata.csv`). Removed references to the nonexistent `interim/mash_scrubbed_species_metadata_2b.csv`.
3. **Cell 3 (MLST parsing)** — Replaced fragile `.str.split('/', expand=True)[3]` with `os.path.basename(x).replace('.fna', '')` to extract genome IDs regardless of path depth.
4. **Cell 8 (output)** — Changed output path from `WORKDIR + 'interim/enriched_metadata_2d.csv'` to `os.path.join(temp_folder, '2d_enriched_metadata.csv')`.

## Access Instructions

- **Notebook**: `examples/2d_enrich_metadata.ipynb`
- **Input**: `temp/2b_genome_metadata.csv` (from notebook 2b), `data/processed/mlst_report.txt` (from Snakemake MLST workflow)
- **Output**: `temp/2d_enriched_metadata.csv`
- **Run**: Execute after notebooks 1a–2c and the anno/MLST/CD-HIT Snakemake workflow have completed
