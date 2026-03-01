# Progress Log 11: CAR & Heaps Notebooks

## Align notebooks 3a/3b with 2c/2d directory standards
**Date:** 2026-02-23 22:24 PST

### Summary
Updated `3a_extract_CAR_genomes.ipynb` and `3b_heaps_plot.ipynb` to use the same
`temp_folder`/`data_dir` path conventions established in notebooks 2c and 2d.

### Steps performed

1. **3a_extract_CAR_genomes.ipynb**
   - Cleaned imports: removed `pickle`, `scipy`, `tqdm`, `plotly`, duplicate `pandas`, commented-out ML imports
   - Added plotting setup (`dpi=200`, seaborn palette/context/style) to imports cell
   - Added `temp_folder` and `data_dir` to config cell
   - P matrix path: `WORKDIR/processed/cd-hit-results/` -> `data_dir/processed/cd-hit-results/`
   - Metadata path: `WORKDIR/interim/mash_scrubbed_species_metadata_2b.csv` -> `temp_folder/2d_enriched_metadata.csv`
   - Save path: `WORKDIR/processed/CAR_genomes/` -> `data_dir/processed/CAR_genomes/`
   - Deleted 5 empty trailing cells (cells 10-14)

2. **3b_heaps_plot.ipynb**
   - Cleaned imports: removed `pickle`, `scipy`, `plotly`
   - Merged plotting setup (from cell 4) and `random.seed(42)` (from cell 3) into imports cell
   - Added `temp_folder` and `data_dir` to config cell
   - Deleted old cells 3 and 4 (now merged)
   - P matrix path: `WORKDIR/processed/cd-hit-results/` -> `data_dir/processed/cd-hit-results/`
   - Metadata path: `WORKDIR/interim/mash_scrubbed_species_metadata_2b.csv` -> `temp_folder/2d_enriched_metadata.csv`

### Path convention (all notebooks 2c-3b)
- `temp_folder` = `CONFIG.get("REUSE_TEMP_DIR", "../temp/")` — intermediate CSVs
- `data_dir` = `CONFIG.get("SNAKEMAKE_DATA_DIR", "data/")` — Snakemake outputs (bakta, cd-hit, mash)

## Fix P-matrix/metadata genome ID mismatch
**Date:** 2026-02-23 22:48 PST

### Problem
3a cell 4 crashed with `KeyError`: 18 genome IDs (all `197.xxxxx`) in the metadata (373 rows) were not columns in the P matrix (355 columns). The BAKTA annotation workflow was run against an older genome list, so 18 of the current 373 QC-passing genomes were never annotated. The P matrix (built in 2c from BAKTA output) only contains the 355 annotated genomes. Notebook 2d saves all 373 metadata rows without filtering to the P matrix.

### Fix
Added a reconciliation line after loading metadata in both 3a (cell 3) and 3b (cell 4):
```python
metadata = metadata[metadata['genome_id'].isin(df_genes.columns)]
```
This filters metadata to only genomes present in the P matrix at analysis time, keeping the saved metadata file as a complete audit record.

### Note
To include all 373 genomes, re-run the BAKTA/CD-HIT Snakemake workflow against the current genome list.
