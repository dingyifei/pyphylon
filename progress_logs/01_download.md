# 01 - Switch Pipeline to Campylobacter jejuni via BV-BRC API

**Date**: 2026-02-16 ~22:00 PST

## Summary

Switched the example pipeline from *Streptococcus pyogenes* (static TSV files) to *Campylobacter jejuni* (live BV-BRC API queries). Added a new `query_bvbrc_genomes()` function to `pyphylon/downloads.py` and updated notebook 1a to use it.

## Key Steps

1. Added `query_bvbrc_genomes()` to `pyphylon/downloads.py` — queries the BV-BRC Data API (`/api/genome/`) using RQL syntax with optional filters for `genome_status` and `genome_quality`, supports pagination via `limit(N,offset)`
2. Updated `examples/config.yml`: `SPECIES_NAME` → "Campylobacter jejuni", `PG_NAME` → "CJejuni", added `TAXON_ID: 197`
3. Updated notebook `1a_filter_genomes_for_download.ipynb`:
   - Added `query_bvbrc_genomes` to imports
   - Replaced TSV file loading (cells 4-5) with API query using `TAXON_ID` from config
   - Updated comments in cells 13, 16, 17 for C. jejuni specifics
4. No changes needed to notebook 1b — it reads CSV outputs from 1a and is organism-agnostic
5. Verified: API returns 594 Complete+Good C. jejuni genomes, all existing pytest tests pass
6. **Bugfix**: Initial implementation used `eq(taxon_id,197)` which only matched genomes tagged with exactly taxon_id 197 (522 genomes), missing 72 genomes assigned to subspecies/strain-level taxon IDs (e.g., *C. jejuni subsp. jejuni* = 32022, *C. jejuni subsp. doylei* = 32021, plus 27 strain-specific IDs). Fixed by switching to `eq(taxon_lineage_ids,197)` which queries by lineage and includes all descendant taxa.

## Access

### Run notebook 1a
```bash
# In the examples/ directory (or Docker container)
jupyter notebook 1a_filter_genomes_for_download.ipynb
```

### Test the API function directly
```python
from pyphylon.downloads import query_bvbrc_genomes
df = query_bvbrc_genomes(197, genome_status='Complete', genome_quality='Good')
print(df.shape)  # ~594 genomes
```

### Run existing tests
```bash
pytest pyphylon/test/test_downloads.py -v
```

---

# 01b - Replace Selenium N50 Scraping with NCBI Datasets API

**Date**: 2026-02-16 ~22:25 PST

## Summary

Replaced `get_scaffold_n50_for_species()` in `pyphylon/downloads.py` which used Selenium + Chrome to scrape NCBI for reference genome Scaffold N50. The old approach failed in WSL (chromedriver exits with status 127 — no Chrome/Chromium installed). Replaced with a single `requests.get()` call to the NCBI Datasets v2 REST API.

## Key Steps

1. Rewrote `get_scaffold_n50_for_species()` to use NCBI Datasets v2 API: `GET https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon_id}/dataset_report?filters.reference_only=true&page_size=1` — response path: `reports[0].assembly_stats.scaffold_n50`
2. Deleted 3 Selenium helper functions that were only used by the old implementation: `get_reference_genome_link()`, `get_scaffold_n50()`, `_convert_to_int()`
3. Removed Selenium/BeautifulSoup imports (`selenium`, `webdriver_manager`, `bs4`) from `downloads.py`. Left `selenium`, `webdriver-manager`, `beautifulsoup4` in `requirements.txt`/`conda/environment.yml` (separate cleanup decision)
4. Updated test fixture in `test_downloads.py`: E. coli K12 expected value changed from `4600000` (rounded from "4.6 Mb" string) to `4641652` (exact API value in bp)
5. Verified: all 4 download tests pass, API returns C. jejuni N50 = 1641481, E. coli N50 = 4641652

## Access

### Test the new function
```python
from pyphylon.downloads import get_scaffold_n50_for_species
get_scaffold_n50_for_species(197)   # C. jejuni → 1641481
get_scaffold_n50_for_species(562)   # E. coli  → 4641652
```

### Run tests
```bash
pytest pyphylon/test/test_downloads.py -v
```

---

## 01c - Pin Dependency Versions (numba/numpy Fix)

**Date**: 2026-02-17 ~00:00 PST

### Summary

Fixed `ImportError: Numba needs NumPy 2.3 or less. Got NumPy 2.4.` which prevented `test_models.py` from being collected by pytest. Pinned all dependency versions in `requirements.txt` to the working conda environment.

### Key Steps

1. Diagnosed: `umap-learn` → `numba 0.63.1` requires numpy <2.4, but unpinned `requirements.txt` allowed numpy 2.4.2
2. Downgraded numpy in WSL conda env: `conda install -n pangenome 'numpy<2.4' -y` → installed numpy 2.3.5
3. Pinned all 18 packages in `requirements.txt` to exact versions from the working env (e.g., `numpy==2.3.5`, `pandas==3.0.0`, `scipy==1.17.0`, etc.)
4. Verified: all 22 pytest tests pass (including `test_models.py` which imports umap/numba)

### Access

```bash
# Verify tests pass
MSYS_NO_PATHCONV=1 wsl -d pangenome -- bash -lc "cd /mnt/f/lab_projects/pangenomics/pyphylon && conda run -n pangenome pytest -v"
```

---

## 01d - Configurable CheckM Filtering by Genome Status

**Date**: 2026-02-17 ~01:00 PST

### Summary

Fixed a crash in `filter_by_genome_quality()` when all genomes are "Complete" (no WGS). CheckM filtering was hardcoded to WGS-only, so it passed an empty DataFrame to `_get_kneebow_cutoff()` which raised a `ValueError`. Added a `checkm_filter_statuses` parameter to control which genome statuses get CheckM-filtered, plus empty-DataFrame guards in all helper functions.

### Key Steps

1. Added `checkm_filter_statuses=('WGS',)` parameter to `filter_by_genome_quality()` in `pyphylon/qcqa.py` — default preserves backward compatibility
2. Restructured function body: after L50/N50 and contig filtering, merges Complete + WGS, then splits into `checkm_subset` (statuses matching the parameter) and `skip_subset`, applies CheckM only to the subset, then recombines
3. Added `if .empty: return` early-exit guards to `_filter_by_contig()`, `_filter_checkM_contamination()`, `_filter_checkM_completeness()`, and `_get_kneebow_cutoff()` (defense-in-depth)
4. Updated notebook `1a_filter_genomes_for_download.ipynb` cell 14: passes `checkm_filter_statuses=('Complete', 'WGS')` so CheckM filtering applies to all genome types
5. Verified: both existing qcqa tests pass unchanged (`test_filter_by_species`, `test_filter_by_genome_quality`)

### Access

```python
# Default (backward compatible) — CheckM only on WGS
qcqa.filter_by_genome_quality(species_summary, min_thresh_n50=min_thresh_n50)

# Filter all genome types by CheckM (use for all-Complete datasets)
qcqa.filter_by_genome_quality(
    species_summary,
    min_thresh_n50=min_thresh_n50,
    checkm_filter_statuses=('Complete', 'WGS'),
)

# Apply CheckM to everything (no status filter)
qcqa.filter_by_genome_quality(
    species_summary,
    min_thresh_n50=min_thresh_n50,
    checkm_filter_statuses=None,
)
```

```bash
# Run tests
MSYS_NO_PATHCONV=1 wsl -d pangenome -- bash -lc "cd /mnt/f/lab_projects/pangenomics/pyphylon && conda activate pangenome && pytest pyphylon/test/test_qcqa.py -v"
```

---

## 01e - Fix CheckM Filter Eliminating All Genomes

**Date**: 2026-02-17 ~02:07 PST

### Summary

Fixed CheckM contamination/completeness filter eliminating all 578 genomes that survived L50/N50 filtering. Root cause: kneebow auto-threshold is too aggressive for uniform Complete+Good datasets from BV-BRC (contamination scores clustered at 1.8–2.8%), causing the elbow detection to pick a near-zero cutoff. Also fixed a `contig_n50` column assignment bug.

### Key Steps

1. **Root cause**: `contamination_cutoff=None` + `checkm_filter_statuses=('Complete', 'WGS')` sent all 578 Complete genomes through kneebow auto-detection. Since BV-BRC Complete+Good genomes have uniformly low contamination, kneebow found an extremely low elbow, and the strict `<` comparison filtered everything out.
2. **Fix 1 — Explicit CheckM cutoffs** (`examples/1a_filter_genomes_for_download.ipynb` cell 14): Changed from `contamination_cutoff=None, completeness_cutoff=None` to `contamination_cutoff=5.0, completeness_cutoff=90.0`. These are standard thresholds from Parks et al. (2015) — *CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes* (Genome Research, 25:1043-1055, DOI: 10.1101/gr.186072.114), which defines "near-complete" genomes as >90% completeness and <5% contamination.
3. **Fix 2 — `contig_n50` column bug** (`pyphylon/qcqa.py` line 114): `filtered_species_summary['contig_n50'] = filtered_species_summary['contig_l50'].astype('int')` was assigning `contig_l50` values to `contig_n50`. Fixed to use `filtered_species_summary['contig_n50'].astype('int')`.
4. Verified: all 2 qcqa tests pass (`test_filter_by_species`, `test_filter_by_genome_quality`)

### Access

```bash
# Run tests
MSYS_NO_PATHCONV=1 wsl -d pangenome -- bash -lc "cd /mnt/f/lab_projects/pangenomics/pyphylon && conda run -n pangenome pytest pyphylon/test/test_qcqa.py -v"
```

---

## 01f - Investigation: Anomalous Pattern in Unfiltered Strain Plot

**Date**: 2026-02-17 ~02:30 PST

### Summary

Investigated an anomalous pattern in the unfiltered strain plot (genome_length vs. patric_cds) from notebook 1a. A dense vertical column at genome_length ≈ 1.671 Mb and a tight cluster at ≈ 1.729 Mb sit above the main trend line due to contamination-inflated gene counts. All anomalous genomes trace to a single BioProject and are correctly removed by CheckM filtering.

### Key Steps

1. **Identified anomaly**: In the unfiltered strain plot, 63 genomes form two distinct clusters well above the main genome_length-vs-CDS trend line, suggesting inflated gene counts relative to genome size
2. **Source**: All 63 genomes belong to BioProject PRJNA942088 — a One Health Surveillance study of *C. jejuni* in the Palestinian poultry supply chain (Swiss Tropical and Public Health Institute)
3. **Two batches with near-identical genome lengths**:
   - **PS_2021** (47 genomes): genome_length std = 43.5 bp, contamination 10.6–21.9%
   - **PS_2022** (16 genomes): genome_length std = 7.3 bp, contamination 10.5–13.5%
4. **Why anomalous**: The extremely low genome length standard deviations (43.5 bp and 7.3 bp across dozens of genomes) suggest these are near-identical backbone assemblies rather than independent assemblies. High contamination (10.5–21.9%) inflates the CDS count, pushing them above the main trend line
5. **Filtering outcome**: All 63 genomes removed by the CheckM contamination filter (>5% contamination cutoff), accounting for 51.6% of all 122 removed genomes
6. **Conclusion**: The filtering pipeline is working correctly — no code changes needed. The anomalous pattern in the unfiltered plot is an expected artifact of contaminated assemblies from a single large submission

### Access

```python
# Reproduce the investigation
import pandas as pd
df = pd.read_csv('output/genome_metadata_1a.csv')
prjna = df[df['bioproject_accession'] == 'PRJNA942088']
print(f"PRJNA942088 genomes: {len(prjna)}")
print(f"Contamination range: {prjna['checkm_contamination'].min():.1f}–{prjna['checkm_contamination'].max():.1f}%")

# Check the two batches
for prefix in ['PS_2021', 'PS_2022']:
    batch = prjna[prjna['genome_name'].str.contains(prefix)]
    print(f"{prefix}: n={len(batch)}, genome_length std={batch['genome_length'].std():.1f} bp")
```

---

## 01g - Create cleanup.sh Helper Script

**Date**: 2026-02-17 ~02:43 PST

### Summary

Created `cleanup.sh` in the project root to selectively remove intermediate files from `temp/` and `output/` directories. Supports phase-based filtering (e.g., `-p 1` for all phase-1 files), directory targeting (`-d output`), protected-file skipping, dry-run previews, and force mode.

### Key Steps

1. Created `cleanup.sh` with argument parsing for `-p <phase>`, `-d <dir>`, `--force`, `--dry-run`, and `-h/--help`
2. Phase matching: single digit (`-p 1`) matches `1[a-zA-Z]_*`; digit+letter (`-p 1a`) matches `1a_*`
3. Protected files: items with "protected" in their name (case-insensitive) are skipped by default, only removed with `--force`
4. `readme.md` is always preserved regardless of flags
5. Uses `find -maxdepth 1 -mindepth 1` for both files and subdirectories; `rm -rf` for directories, `rm -f` for files
6. Verified all 5 test cases: full dry-run, phase-1 dry-run, protected-file skip, force mode, temp-only targeting

### Access

```bash
# Preview all deletions
./cleanup.sh --dry-run

# Clean all phase-1 files (1a_*, 1b_*, ...)
./cleanup.sh -p 1

# Clean only 1a_* in output/
./cleanup.sh -p 1a -d output

# Include protected files
./cleanup.sh --force

# Preview phase-2 cleanup
./cleanup.sh -p 2 --dry-run

# Show help
./cleanup.sh -h
```

---

## Notes

- BV-BRC API endpoint: `https://www.bv-brc.org/api/genome/`
- Query language: RQL (Resource Query Language) — `eq(field,value)&limit(N,offset)`
- Response format: JSON array with `Accept: application/json` header
- Returns ~99 columns per genome including genome_id, genome_name, genome_length, patric_cds, contigs, contig_n50, gc_content, checkm_completeness, checkm_contamination, etc.
- Column names match the BV-BRC TSV files, so downstream processing is compatible
- The number of genomes may change over time as BV-BRC is updated (594 as of 2026-02-17)
- **Important**: Use `eq(taxon_lineage_ids,ID)` not `eq(taxon_id,ID)` when querying by taxonomy. The latter only matches the exact taxon ID and misses subspecies/strain-level assignments. `taxon_lineage_ids` searches the full lineage and includes all descendants.
