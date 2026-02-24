# 02: Fix CheckM filtering dropping genomes with null CheckM data

**Date**: 2026-02-17, 2:24 AM PST

## Problem

`filter_by_genome_quality()` silently dropped genomes with null CheckM data via `dropna()` inside the helper functions. For *C. jejuni*, 22 Complete+Good genomes (a clean Sept-2022 BV-BRC batch) had no CheckM scores and were lost — going from 578 post-contig-filter to 471 instead of ~493.

## Key Steps

1. **Added `checkm_missing` parameter** to `filter_by_genome_quality()` in `pyphylon/qcqa.py`
   - `'keep'` (new default) — retains genomes with null CheckM values
   - `'drop'` — discards them (old behavior)

2. **Rewrote `_filter_checkM_contamination()` and `_filter_checkM_completeness()`** to split data into has_data/missing_data, filter only has_data, then recombine based on `checkm_missing` setting

3. **Fixed post-filter typecast** — changed `astype('float')` to `pd.to_numeric(..., errors='coerce')` for `checkm_contamination` and `checkm_completeness` columns to handle NaN values

4. **Updated notebook 1a cell 14** — added `checkm_missing='keep'` to the call

5. **Updated test** in `pyphylon/test/test_qcqa.py` — tests both `'keep'` (1808 genomes) and `'drop'` (257 genomes) modes

6. **Tests pass**: `pytest pyphylon/test/test_qcqa.py -v` — 2/2 passed

## Files Modified

- `pyphylon/qcqa.py` — new parameter, rewritten CheckM helpers, safe typecasts
- `pyphylon/test/test_qcqa.py` — updated assertions for both modes
- `examples/1a_filter_genomes_for_download.ipynb` — cell 14 updated

## Expected Result

Re-running notebook 1a should yield ~493 genomes (471 passing CheckM + 22 null retained) instead of 471.

---

## Re-run: Feb 23, 2026

Notebook 1a was re-run with two parameter changes from the original run.

### Parameter Changes

| Parameter | Previous | Re-run |
|-----------|----------|--------|
| `checkm_missing` | `'keep'` | `'drop'` — now drops genomes with null CheckM data |
| Completeness plot cutoff line | `y=90` | `y=92` — corrected to match actual `completeness_cutoff=92.0` |

### Filtration Results

| Stage | Previous run | Re-run |
|-------|:---:|:---:|
| BV-BRC query (Complete+Good) | 594 | 594 |
| Species filter | 593 | 593 |
| L50/N50 | 578 | 578 |
| CheckM filtered out | 85 | 109 |
| **After notebook 1a** | **493** | **469** |

The 24-genome difference (493 → 469) comes from the 22 null-CheckM genomes now being dropped plus 2 additional genomes caught by the stricter completeness cutoff (92% vs 90%).

### Note

Downstream notebooks (1b, 2a, 2b, etc.) have not been re-run yet with the updated 469-genome set.
