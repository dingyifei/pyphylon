# 06 — Download PDFs for Papers in papers.md

**Date**: 2026-02-17 04:19 PST

## Summary

Automated download of open-access PDFs for the 55 unique papers listed in `papers.md`. Retrieved 46/55 PDFs (84 MB) and annotated `papers.md` with `[PDF]` links.

## Key Steps

1. **Created `papers/` directory** and added `papers/*.pdf` to `.gitignore`
2. **Wrote `download_papers.py`** — self-contained download script with:
   - Regex parsing of `papers.md` to extract DOI, PMID, arXiv ID, author, year
   - Tiered URL resolution:
     - Tier 1: Europe PMC direct PDF (`europepmc.org`) for papers with PMCIDs (42 found via NCBI ID Converter)
     - Tier 2: Unpaywall API (`api.unpaywall.org`) for publisher OA / repository copies
     - Tier 3: arXiv direct for arXiv-only entries
   - Parallel downloads (ThreadPoolExecutor, 8 workers) with `%PDF` magic byte validation
   - Automatic `papers.md` update — appends `| [PDF](papers/Author_Year.pdf)` to each retrieved entry
   - CLI flags: `--verbose`, `--dry-run`, `--retry-failed`
3. **Debugged and iterated** through 3 issues:
   - Regex: initial single-regex approach failed; switched to independent DOI/PMID/arXiv extraction per line
   - NCBI ID Converter: response returns PMID as int, not string — fixed with `str()` cast
   - PMC 403 errors: NCBI PMC now blocks programmatic PDF access; switched to Europe PMC endpoint
4. **Ran final download** — 46 PDFs retrieved, 50 `[PDF]` links added to `papers.md` (46 unique + 4 §8 duplicates)

## Results

| Category | Count | Notes |
|----------|-------|-------|
| Downloaded | 46 | Valid PDFs, `%PDF` magic byte verified |
| Paywalled | 5 | Elsevier (3), Nature Reviews (1), OUP (1) — no OA copies |
| Failed | 4 | Wiley (2), Nature (1), Wiley/repository (1) — URL resolved but PDF blocked |
| Duplicates (§8) | 4 | Cross-references linked to same PDFs as §1-7 |

### Paywalled (no OA available)

- Sheppard et al. (2018) — Elsevier book chapter
- Young et al. (2007) — Nature Reviews Microbiology
- Dasti et al. (2010) — Elsevier IJMM
- Bolton (2015) — Elsevier Food Microbiology
- Seemann (2014) — OUP Bioinformatics

### Failed (URL resolved, download blocked)

- Yahara et al. (2017) — repository landing page, no direct PDF
- Baily et al. (2015) — repository landing page, no direct PDF
- Hepworth et al. (2011) — Wiley paywall
- Burnham & Hendrixson (2018) — Nature paywall

## Access Instructions

```bash
# Run the download script
python download_papers.py --verbose

# Re-attempt previously failed downloads
python download_papers.py --retry-failed --verbose

# Preview what would be downloaded (no actual downloads)
python download_papers.py --dry-run

# Check download report
cat papers/_download_report.json

# Verify all PDFs are valid
file papers/*.pdf | grep -v "PDF document"
```

## Files Created/Modified

| File | Action |
|------|--------|
| `download_papers.py` | Created — download script (~300 lines) |
| `papers/` | Created — 46 PDFs (84 MB total) |
| `papers/_download_report.json` | Created — full download results |
| `papers.md` | Modified — 50 `[PDF]` links added |
| `.gitignore` | Modified — added `papers/*.pdf` |
