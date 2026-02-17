# 03 — Fix BV-BRC Genome Downloads (FTP → HTTPS)

**Date**: 2026-02-17 03:15 PST

## Problem

ALL genome downloads in notebook 1b failed with "Bad genome ID" warnings, resulting in 0 genomes downloaded. Root cause: BV-BRC's FTP server now enforces SSL/TLS on the control channel, and Python's `urllib.request.urlretrieve()` uses plain FTP, which gets rejected with `ftplib.error_perm: 550 SSL/TLS required on the control channel`. The `IOError` catch misleadingly logged every failure as "Bad genome ID".

## Changes Made

### 1. `pyphylon/downloads.py` — Switch FTP to HTTPS Data API

- **Added `BV_BRC_API_BASE` and `BV_BRC_API_ENDPOINTS` constants** (after `VALID_BV_BRC_FILES`): Maps each filetype (`fna`, `gff`, `faa`, `ffn`, `frn`, `features.tab`, `pathway.tab`, `spgene.tab`, `subsystem.tab`) to its BV-BRC API resource and `Accept` header
- **Replaced FTP download block in `download_genomes_bvbrc()`**: Swapped `urllib.request.urlretrieve()` (FTP) for `urllib.request.urlopen()` (HTTPS) using the BV-BRC Data API endpoints. URL pattern: `GET /api/{resource}/?eq(genome_id,{id})&http_accept={header}&limit(25000)`
- **Improved error handling**: Catches `urllib.error.HTTPError` / `URLError` instead of bare `IOError`; logs descriptive messages (`"Failed to download {filetype} for genome {id}: {error}"`) instead of `"Bad genome ID"`
- **Added empty-response guard**: Warns if API returns empty data and marks genome as bad
- **Fixed cleanup bug (line 313)**: Changed `f"{genome}.{ftype}"` → `f"{bad_genome}.{ftype}"` — the old code used a stale loop variable `genome` from the download loop, targeting the wrong files for deletion

### 2. `examples/1b_download_genomes_bvbrc.ipynb` — Notebook modernization

- **Added `yaml` import** and `config.yml` loading for `REUSE_TEMP`, `REUSE_TEMP_DIR`, `DEBUG` settings
- **Changed download directory**: `data/raw/genomes/` → `temp/1b_protected/raw/genomes/` (survives `cleanup.sh` unless `--force`)
- **Updated input paths**: Reads from `output/1a_genome_summary.csv` and `output/1a_genome_metadata.csv` (output of notebook 1a) instead of hardcoded `data/interim/` paths
- **Updated output paths**: Saves to `output/1b_genome_summary.csv` and `output/1b_genome_metadata.csv` instead of `data/interim/`
- **Removed empty trailing cells**

### 3. `examples/config.yml` — New config file (untracked)

- Added `REUSE_TEMP: True` and `REUSE_TEMP_DIR: "../temp/"` settings
- Added `DEBUG: True` flag
- Contains `METADATA_FILE`, `GENOMES_FILE`, `FILTERED_GENOMES_FILE` path settings

### 4. `.gitignore` — Minor updates

- Added `/.wsl` (WSL disk images) and `/.idea` (JetBrains IDE)
- Minor formatting fix (lost `# Byte-` comment prefix — cosmetic)

### 5. `input/readme.md` and `temp/readme.md` — Added placeholder content

- `input/readme.md`: "Put input files here"
- `temp/readme.md`: Describes temp folder purpose and `reuse_temp` flag

## Verification

- All 4 tests in `test_downloads.py` pass
- `BV_BRC_API_BASE` and `BV_BRC_API_ENDPOINTS` import correctly
- HTTPS API verified to return valid FASTA for C. jejuni genomes (~1.8 MB each)

## Key API Details

| Filetype | API Resource | Accept Header |
|----------|-------------|---------------|
| `fna` | `genome_sequence` | `application/sralign+dna+fasta` |
| `gff` | `genome_feature` | `application/gff` |
| `faa` | `genome_feature` | `application/protein+fasta` |
| `ffn` | `genome_feature` | `application/dna+fasta` |
| `frn` | `genome_feature` | `application/dna+fasta` |
| `features.tab` | `genome_feature` | `text/tsv` |
| `pathway.tab` | `pathway` | `text/tsv` |
| `spgene.tab` | `sp_gene` | `text/tsv` |
| `subsystem.tab` | `subsystem` | `text/tsv` |