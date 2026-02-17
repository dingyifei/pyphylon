# 04 - Update 2a_clean_metadata.ipynb

**Date**: 2026-02-17 03:56 AM PST

## Summary

Updated notebook `examples/2a_clean_metadata.ipynb` to work with the updated 1a/1b pipeline outputs and fixed a summary sync bug.

## Key Steps

1. **Updated imports and config loading** (cell 2) — Added `os`, `load_config` imports and `output_folder` path variable, matching the 1b notebook pattern
2. **Updated input paths** (cell 3) — Changed reads from `data/interim/genome_*_1b.csv` to `../output/1b_genome_*.csv` (where 1b actually writes)
3. **Fixed summary sync bug** (cell 6) — After deduplicating metadata by `biosample_accession`, added `downloaded_species_summary = downloaded_species_summary.loc[downloaded_species_metadata.index]` to keep summary in sync
4. **Updated output paths** (cell 12) — Changed writes from `data/interim/genome_*_2a.csv` to `../output/2a_genome_*.csv` (following `{step}_{type}.csv` convention)

## Files Modified

- `examples/2a_clean_metadata.ipynb` — 4 cells updated (cells 2, 3, 6, 12)

## Access Instructions

- Run the notebook from the `examples/` directory: `jupyter notebook examples/2a_clean_metadata.ipynb`
- Requires `output/1b_genome_summary.csv` and `output/1b_genome_metadata.csv` from notebook 1b
- Outputs: `output/2a_genome_summary.csv` and `output/2a_genome_metadata.csv`

---

# 04b - Update 2b_mash_filtration_and_clustering.ipynb

**Date**: 2026-02-17 04:02 AM PST

## Summary

Updated notebook `examples/2b_mash_filtration_and_clustering.ipynb` to read 2a outputs from `../output/` and write its own outputs there too. Snakemake workflow paths (raw genomes, `mash_distances.txt`) remain under `WORKDIR`.

## Key Steps

1. **Added `output_folder`** (cell 2) — Added `output_folder = os.path.join("../output/")` alongside existing `WORKDIR`
2. **Updated input paths** (cell 4) — Changed reads from `WORKDIR/interim/genome_*_2a.csv` to `output_folder/2a_genome_*.csv`
3. **Updated correlation matrix cache path** (cell 11) — Changed from `WORKDIR + 'processed/df_mash_corr_dist.csv'` to `output_folder/2b_mash_corr_dist.csv`
4. **Updated summary output** (cell 43) — Changed from `WORKDIR/interim/mash_scrubbed_species_summary_2b.csv` to `output_folder/2b_genome_summary.csv`
5. **Updated metadata output** (cell 44) — Changed from `WORKDIR/interim/mash_scrubbed_species_metadata_2b.csv` to `output_folder/2b_genome_metadata.csv`
6. **Updated mash matrix outputs** (cell 46) — Changed from `WORKDIR/interim/df_mash_*.csv` to `output_folder/2b_mash_square.csv` and `output_folder/2b_mash_corr_dist.csv`
7. **Updated config.yml** — Changed `FILTERED_GENOMES_FILE` from `data/interim/mash_scrubbed_species_metadata_2b.csv` to `output/2b_genome_metadata.csv`

## Files Modified

- `examples/2b_mash_filtration_and_clustering.ipynb` — 6 cells updated (cells 2, 4, 11, 43, 44, 46)
- `examples/config.yml` — `FILTERED_GENOMES_FILE` path updated

## Access Instructions

- Run the notebook from the `examples/` directory: `jupyter notebook examples/2b_mash_filtration_and_clustering.ipynb`
- Requires `output/2a_genome_summary.csv` and `output/2a_genome_metadata.csv` from notebook 2a
- Requires Mash workflow to have been run first (reads `WORKDIR/processed/mash/mash_distances.txt`)
- Outputs: `output/2b_genome_summary.csv`, `output/2b_genome_metadata.csv`, `output/2b_mash_square.csv`, `output/2b_mash_corr_dist.csv`

## Downstream Note

Notebooks 2c and 3a still hardcode `WORKDIR + 'interim/mash_scrubbed_species_metadata_2b.csv'` and will need separate updates to read from `../output/2b_genome_metadata.csv`.
