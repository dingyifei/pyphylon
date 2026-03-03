# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Pyphylon is a Python bioinformatics package for analyzing and visualizing co-occurring gene/allele sets (phylons) within pangenomes. It combines automated genomic data workflows with advanced computational analysis methods including Non-negative Matrix Factorization (NMF), Multiple Correspondence Analysis (MCA), UMAP, and HDBSCAN clustering.

## Development Commands

### Pipeline Modes

The pipeline supports two execution modes. The Snakefile's `container:` directives are ignored when `--sdm apptainer` is not passed, so the same workflow works in both modes.

**Conda native** (macOS ARM / Linux) — all bio tools installed via conda, no Docker needed:
```bash
conda env update -f conda/environment-marimo.yml
conda activate pyphylon-marimo
snakemake --cores 8 --resources jupyter_kernel=1
```

**Docker** (amd64 platforms / CI) — containerized bio tools via Apptainer:
```bash
docker compose build
docker compose run --rm pipeline snakemake --cores 8 --sdm apptainer
```

### Environment Setup (macOS)

**Conda environment**: `pyphylon-marimo` (miniconda, Python 3.x, marimo 0.20.2, pandas 3.0.0)

```bash
# Activate the environment
conda activate pyphylon-marimo

# Run a marimo notebook headlessly
conda run -n pyphylon-marimo python notebooks/1a_filter_genomes.py -- --config config.yml

# Run a marimo notebook interactively
conda run -n pyphylon-marimo marimo edit notebooks/1a_filter_genomes.py

# Lint and format
conda run -n pyphylon-marimo ruff check notebooks/
conda run -n pyphylon-marimo ruff format notebooks/
```

### Cleanup
```bash
# Remove all intermediate files from temp/ and output/ (preserves readme.md + protected items)
./cleanup.sh

# Preview what would be deleted
./cleanup.sh --dry-run

# Clean specific phase (e.g., all phase-1 files: 1a_*, 1b_*, ...)
./cleanup.sh -p 1

# Clean specific phase+step in one directory
./cleanup.sh -p 1a -d output

# Include files with "protected" in their name
./cleanup.sh --force
```

### Testing
```bash
# Run all tests
pytest

# Run specific test modules
pytest pyphylon/test/test_models.py
pytest pyphylon/test/test_pangenome.py
pytest pyphylon/test/test_downloads.py
pytest pyphylon/test/test_qcqa.py
pytest pyphylon/test/test_util.py
```

### Docker Workflows

Uses `docker-compose.yml` with OrbStack (macOS) or Docker Desktop. The entire repo is bind-mounted; code changes are instantly visible without rebuilding.

```bash
# Build image (only needed when requirements.txt changes)
docker compose build

# Full pipeline (notebooks + bio workflows)
docker compose run pipeline snakemake --cores 4 --sdm apptainer

# Bio workflow only
docker compose run pipeline snakemake pangenome_mash --cores 4 --sdm apptainer

# Single notebook rule
docker compose run pipeline snakemake nb_1a --cores 4

# Interactive shell
docker compose run pipeline bash
```

## Architecture Overview

### Core Design Pattern
The package follows a three-tier architecture centered around the `NmfData` object (`pyphylon/core.py`), which serves as the primary analysis container:

1. **Data Acquisition & Preprocessing Layer**
   - `downloads.py`: Automated genome downloads from BV-BRC and NCBI
   - `qcqa.py`: Quality control and assurance functions
   - Snakemake workflows for BAKTA annotation, MLST typing, Mash clustering

2. **Pangenome Analysis Layer**
   - `pangenome.py`: Core pangenome construction using CD-HIT clustering
   - `mash.py`: Genome distance metrics and sketching
   - Matrix building and validation (P, L, A, V, U, F matrices)

3. **Computational Analysis & Visualization Layer**
   - `models.py`: NMF, MCA, UMAP, HDBSCAN implementations
   - `plotting.py`: Circular genome plots, dendrograms, visualizations
   - `biointerp.py`: Biological interpretation tools
   - `io.py`: Data serialization using joblib

### Key Data Structures
The `NmfData` class in `core.py` manages all analysis matrices:
- **P matrix**: Presence/absence matrix (strains × genes)
- **L matrix**: Left factorization matrix (phylons × strains)
- **A matrix**: Right factorization matrix (phylons × genes)
- **V, U, F matrices**: Various normalized/binary transformations
- **metadata_table**: Strain metadata (MLST, serotype, etc.)

### Workflow Integration
The package integrates marimo notebooks (`notebooks/`) with Snakemake workflows (`workflow/`) for reproducible bioinformatics pipelines. The master `Snakefile` uses `include:` to pull in three bio workflow Snakefiles:

- **Notebooks 1a-1b**: Genome filtering and downloading from BV-BRC
- **Mash workflow** (`workflow/mash/`): Genome distance calculation — feeds into nb_2b
- **Notebooks 2a-2b**: Mash-based strain filtering (2b is a Snakemake checkpoint)
- **Anno/MLST/CD-HIT workflow** (`workflow/anno_mlst_cdhit/`): Annotation, typing, pangenome construction — feeds into nb_2c, nb_2d
- **Notebooks 2c-2d**: Pangenome matrix generation and metadata enrichment
- **Notebooks 3a-4a**: Core/accessory/rare analysis and NMF decomposition
- **Notebooks 5a-5e**: Phylon characterization and functional analysis
- **Infer affinities workflow** (`workflow/infer_affinities/`): Compare new strains to pangenome — feeds into nb_5f
- **Notebook 5f**: Affinity inference visualization

## Key Configuration Files

### Workflow Configuration (`config.yml` — repository root)
```yaml
SPECIES_NAME: "Campylobacter jejuni"
PG_NAME: "CJejuni"
TAXON_ID: 197
DATA_DIR: "data/"
TEMP_DIR: "temp/"
OUTPUT_DIR: "output/"
BAKTA_DB_DIR: "temp/db-light"
BAKTA_THREADS: 8
CDHIT_THREADS: 8
MLST_THREADS: 4
CONTAINER_BAKTA: "docker://oschwengers/bakta@sha256:..."
CONTAINER_MASH: "docker://staphb/mashtree@sha256:..."
CONTAINER_CDHIT: "docker://biocontainers/cd-hit@sha256:..."
CONTAINER_MLST: "docker://staphb/mlst@sha256:..."
```

### BV-BRC API
The pipeline can fetch genome metadata directly from the BV-BRC Data API using `query_bvbrc_genomes()` in `pyphylon/downloads.py`:
- **Endpoint**: `https://www.bv-brc.org/api/genome/`
- **Query format**: RQL — `eq(taxon_id,197)&eq(genome_status,Complete)&eq(genome_quality,Good)&limit(25000,0)`
- **Response**: JSON array → `pd.DataFrame` (same column names as BV-BRC TSV downloads)
- **Usage**: `query_bvbrc_genomes(taxon_id, genome_status=None, genome_quality=None, limit=25000)`

### Dependencies
- **Core**: Python >=3.11, numpy, pandas, scipy, scikit-learn, biopython
- **Analysis**: prince (MCA), hdbscan, umap-learn, kneebow
- **Visualization**: matplotlib, seaborn, plotly
- **Web scraping**: selenium, webdriver-manager, beautifulsoup4
- **Bioinformatics tools**: BAKTA (annotation), MLST, CD-HIT, Mash

## Data Flow and File Structure

### Directory Structure
```
data/                          # Git submodule (LFS-backed)
├── raw/genomes/fna/           # Input genome files (.fna)
├── processed/
│   ├── bakta/                 # Genome annotations
│   ├── cd-hit-results/        # Pangenome clusters
│   ├── mlst/                  # MLST typing results
│   ├── CAR_genomes/           # Core/accessory/rare analysis
│   └── nmf-outputs/           # NMF decomposition results
├── inferring_affinities/      # New strain data for affinity inference
temp/                          # Intermediate files (gitignored)
├── 1a_*, 1b_*, 2a_*, 2b_*    # Per-step CSV intermediates
├── 2b_mash/                   # Mash sketches and distances
└── db-light/                  # BAKTA database
output/
├── figures/                   # PNG plots per step
└── data/                      # Summary CSVs per step
```

### Critical Matrix Files
- `SPyogenes_strain_by_gene.pickle.gz`: Main P matrix (strains × genes)
- `A.csv`, `L.csv`: NMF factorization matrices
- `A_norm.csv`, `L_norm.csv`, `A_bin.csv`, `L_bin.csv`: Normalized/binary variants

## Development Notes

### Testing Strategy
- Test data located in `pyphylon/test/data/bakta/` with 3 sample genomes
- Tests cover models, pangenome construction, downloads, QC/QA, and utilities
- Use pytest framework with modular test files

### External Tool Dependencies
- **BAKTA**: Genome annotation (containerized)
- **MLST**: Multi-locus sequence typing (containerized)
- **CD-HIT**: Sequence clustering for pangenome construction
- **Mash**: Fast genome distance estimation (containerized via staphb/mashtree)

### Serialization and I/O
- Primary serialization uses joblib for large matrix objects
- CSV exports for interoperability with R/other tools
- Metadata handling through pandas DataFrames

### Performance Considerations
- Large genome datasets require careful memory management
- Docker containerization ensures reproducible bioinformatics environments
- Snakemake provides parallel processing capabilities
- Matrix operations optimized through scipy sparse matrices where applicable

## WSL Environments

### pangenome (Fedora 43)
- **Distro name**: `pangenome`
- **Base**: Fedora 43 (cloned from FedoraLinux-43)
- **Disk location**: `.wsl\pangenome\` (project-local)
- **Default user**: `yifei` (passwordless sudo)
- **Access**: `wsl -d pangenome`
- **Docker**: Docker CE 29.2.1 installed, `yifei` in docker group. Start daemon: `sudo dockerd &` or `sudo systemctl start docker`
- **Conda**: Miniforge at `/home/yifei/miniforge3`, env `pangenome` (Python 3.11) with all pyphylon deps + pyphylon installed in editable mode
- **Activate env**: `conda activate pangenome`
- **Docker image**: `pyphylon:latest` (snakemake/snakemake:v8.30.0 base, Debian Bookworm, 3.96GB)
- **Note**: Use `MSYS_NO_PATHCONV=1` prefix when running wsl commands from Git Bash to avoid path mangling.

## Progress Tracking

All 15 notebook conversions (1a through 5f) are complete. The pipeline is orchestrated by the master `Snakefile` with three included bio workflow Snakefiles.