# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Pyphylon is a Python bioinformatics package for analyzing and visualizing co-occurring gene/allele sets (phylons) within pangenomes. It combines automated genomic data workflows with advanced computational analysis methods including Non-negative Matrix Factorization (NMF), Multiple Correspondence Analysis (MCA), UMAP, and HDBSCAN clustering.

## Development Commands

### Environment Setup

**Note**: Development environment uses PowerShell on Windows with WSL containers for Docker/conda workflows.

```bash
# Clone and install in development mode
git clone https://github.com/SBRG/pyphylon.git
cd pyphylon
pip install -r requirements.txt
pip install -e .

# Optional: Install with extras
pip install -e .[cd-hit,tests]

# Using conda environment (miniforge available in WSL)
conda env create -f conda/environment.yml
conda activate pyphylon
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

**Note**: Base shell is PowerShell on Windows. WSL containers have full Docker and miniforge conda available.

```bash
# Build container
docker build -t pyphylon .

# Run bioinformatics workflows
# PowerShell (Windows):
docker run --privileged -it -v ${PWD}/examples:/examples -v ${PWD}/workflow:/workflow pyphylon

# WSL2/Linux:
docker run --privileged -it -v $(pwd)/examples:/examples -v $(pwd)/workflow:/workflow pyphylon

# Inside container - Mash workflow
cd workflow/mash
snakemake -d /examples/data --use-singularity -c 10

# Inside container - Annotation, MLST, CD-HIT workflow
cd /workflow/anno_mlst_cdhit
snakemake -d /examples/data --use-singularity -c 10
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
The package integrates Jupyter notebooks (`examples/`) with Snakemake workflows (`workflow/`) for reproducible bioinformatics pipelines:

- **Notebooks 1a-1b**: Genome filtering and downloading from BV-BRC
- **Mash workflow**: Genome distance calculation and species filtering
- **Notebooks 2a-2b**: Mash-based strain filtering
- **Anno/MLST/CD-HIT workflow**: Annotation, typing, pangenome construction
- **Notebooks 2c-2d**: Pangenome matrix generation and metadata enrichment
- **Notebooks 3a-4a**: Core/accessory/rare analysis and NMF decomposition
- **Notebooks 5a-5f**: Phylon characterization and functional analysis

## Key Configuration Files

### Workflow Configuration (`examples/config.yml`)
```yaml
WORKDIR: "data/"
SPECIES_NAME: "Campylobacter jejuni"
PG_NAME: "CJejuni"
TAXON_ID: 197
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

### Expected Directory Structure for Examples
```
examples/data/
├── raw/genomes/fna/          # Input genome files (.fna)
├── interim/                  # Intermediate CSV files
├── processed/
│   ├── bakta/               # Genome annotations
│   ├── cd-hit-results/      # Pangenome clusters
│   ├── mlst/                # MLST typing results
│   ├── mash/                # Genome sketches/distances
│   ├── CAR_genomes/         # Core/accessory/rare analysis
│   └── nmf-outputs/         # NMF decomposition results
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

## Progress Logs

Progress logs are stored in `progress_logs/` and track work done on this project.

### Format Requirements
Each log entry must include:
- **Date/time** of the operation
- **Key steps** performed (numbered)
- **Access instructions** for anything set up (commands, paths, URLs)