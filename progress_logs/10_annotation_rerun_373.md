# 10 — Rerun BAKTA/MLST/CD-HIT for 373-Genome Subset

**Date**: 2026-02-23 21:52 PST

## Overview

Re-ran the BAKTA annotation, MLST typing, and CD-HIT clustering pipeline for 373 genomes (reduced from 467 by re-enabling `small_clst_limit=5` in notebook 2b). All 373 are a subset of the previously processed 467, so most BAKTA/MLST outputs already existed.

## Key Steps

1. **Started Docker daemon** in WSL pangenome distro
2. **Copied updated `2b_genome_metadata.csv`** (374 lines = 373 genomes + header) from Windows FS → `/home/yifei/pyphylon_run/examples/temp/`
3. **Deleted stale aggregate outputs** on native Linux FS (required `sudo` — files owned by root from Docker):
   - `raw/CJejuni.faa` (concatenated protein sequences)
   - `processed/cd-hit-results/*` (CD-HIT cluster output)
   - `processed/mlst_report.txt` (combined MLST results)
4. **Removed old Docker container** (`pyphylon_anno`)
5. **Launched Snakemake** in Docker container:
   - 40 total jobs: 18 bakta_annotation + 18 mlst + concat_faas + cdhit + mlst_report + all
   - 18 genomes had missing BAKTA/MLST outputs (likely from incomplete runs previously)
   - **Runtime: ~6 minutes** (vs ~2 hours for full 467-genome run)
6. **Verified outputs** on native Linux FS:
   - `CJejuni`: 2.5 MB, `CJejuni.clstr`: 26 MB
   - `mlst_report.txt`: 373 lines (exact match)
   - `CJejuni.faa`: 230 MB
7. **Copied results** back to Windows FS

## Results on Windows FS

```
examples/data/processed/cd-hit-results/
├── CJejuni          (2.5 MB, Feb 23)    ← NEW: 373-genome CD-HIT output
├── CJejuni.clstr    (26 MB, Feb 23)     ← NEW: 373-genome cluster file
├── CJejuni.cdhit.clstr    (3.4 MB, Feb 17)  ← OLD: from notebook 2c (467 genomes)
├── CJejuni_allele_names.tsv             ← OLD: from notebook 2c (467 genomes)
├── CJejuni_strain_by_allele.pickle.gz   ← OLD: from notebook 2c (467 genomes)
├── CJejuni_strain_by_gene.pickle.gz     ← OLD: from notebook 2c (467 genomes)

examples/data/processed/mlst_report.txt  (373 lines, Feb 23)
```

**Note**: Old Feb 17 files (`CJejuni.cdhit.clstr`, `*_allele_names.tsv`, `*.pickle.gz`) are stale from the 467-genome run. They will be overwritten when notebook 2c runs `build_cds_pangenome()` with the new CD-HIT output.

## Access Instructions

```bash
# Verify CD-HIT output
ls -la examples/data/processed/cd-hit-results/

# Verify MLST report (expect 373 lines)
wc -l examples/data/processed/mlst_report.txt
```

## Status

- [x] Docker daemon started
- [x] Updated genome metadata copied to native Linux FS
- [x] Stale aggregate outputs deleted
- [x] Snakemake rerun complete (40/40 steps)
- [x] Results verified on native Linux FS
- [x] Results copied to Windows FS
- [x] Final verification on Windows FS
