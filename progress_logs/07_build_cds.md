# 07 — Update Notebook 2c for Modernized Pipeline

**Date**: 2026-02-17, 5:53 AM PST

## Summary

Updated `examples/2c_build_cds_pangenome.ipynb` to work with the modernized `temp/` + `output/` directory convention established in notebooks 1a-2b. The notebook previously used stale paths that no longer existed after `WORKDIR` was changed from `"data/"` to `"../temp/"`.

## Key Steps

1. **Added `SNAKEMAKE_DATA_DIR: "data/"` to `examples/config.yml`** — separates the Snakemake output directory (BAKTA, CD-HIT results) from the notebook intermediate directory (`temp/`).

2. **Cleaned up imports (cell 0)** — removed unused `import pickle` and duplicate `import pandas as pd`.

3. **Added `temp_folder` and `data_dir` variables (cell 1)** — follows the same `CONFIG.get()` pattern used in notebook 2b.

4. **Deleted unused `import gzip` cell** — was cell 2, never referenced.

5. **Fixed metadata input path (cell 2)** — changed from `{WORKDIR}/interim/mash_scrubbed_species_metadata_2b.csv` (doesn't exist) to `{temp_folder}/2b_genome_metadata.csv` (actual 2b output).

6. **Fixed BAKTA path (cell 3)** — changed from `{WORKDIR}/processed/bakta/` (= `../temp/processed/bakta/`, wrong) to `{data_dir}/processed/bakta/` (= `data/processed/bakta/`, where Snakemake writes). Added `os.path.isdir()` guard to skip non-directory entries.

7. **Rewrote genome filter (cell 5)** — replaced O(n^2) nested loop with O(n) set-based lookup using `os.path.basename(os.path.dirname(f))`.

8. **Fixed CD-HIT output path (cell 7)** — changed from hardcoded `'data/processed/cd-hit-results/'` to `{data_dir}/processed/cd-hit-results/` with `os.makedirs(exist_ok=True)`.

## Files Modified

- `examples/config.yml` — added `SNAKEMAKE_DATA_DIR` key
- `examples/2c_build_cds_pangenome.ipynb` — 6 cells modified, 1 cell deleted

## Final Cell Layout

| Cell | Content |
|------|---------|
| 0 | Imports (cleaned) |
| 1 | Config + `temp_folder`, `data_dir` |
| 2 | Load `2b_genome_metadata.csv` from `temp/` |
| 3 | Build BAKTA `.faa` paths from `data/` |
| 4 | Sanity check (unchanged) |
| 5 | Set-based filter to matching genomes |
| 6 | `len(real_paths)` (unchanged) |
| 7 | `build_cds_pangenome` with config-driven output |
| 8 | `df_genes.sum()` (unchanged) |
| 9 | Clustermap (unchanged) |

## Verification

- Cells 0-6 can be validated once `temp/2b_genome_metadata.csv` and `data/processed/bakta/` exist
- Cell 7 requires CD-HIT installed (available in WSL `pangenome` conda env)
