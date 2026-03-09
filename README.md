# pyphylon

Python package for analyzing and visualizing co-occurring gene/allele sets (phylons) within a pangenome. The pipeline combines **marimo** interactive notebooks, **Snakemake** orchestration, containerized bioinformatics tools, and **Quarto** PDF reports.

## Quick Start

### Baremetal (conda)

```bash
pip install -e .
conda run -n pyphylon-marimo python notebooks/1a_filter_genomes.py -- --config config.yml
```

### Docker (OrbStack)

```bash
docker compose build
docker compose run pipeline snakemake --cores 4 --sdm apptainer
```

## Directory Structure

```
pyphylon-marimo/
├── notebooks/          # Marimo .py notebooks (1a through 5f)
├── reports/            # Quarto .qmd reports
├── workflow/           # Bio Snakefiles (mash, anno_mlst_cdhit, infer_affinities)
├── pyphylon/           # Python package
├── data/               # Git submodule (LFS-backed)
├── temp/               # Intermediate files (gitignored)
├── output/             # Figures, data CSVs, PDF reports
├── Snakefile           # Master orchestration
├── config.yml          # Unified config
└── docker-compose.yml  # Docker pipeline
```

## Running the Pipeline

### Host (development)

```bash
# Run a single notebook headlessly
conda run -n pyphylon-marimo python notebooks/1a_filter_genomes.py -- --config config.yml

# Run a notebook interactively in the browser
conda run -n pyphylon-marimo marimo edit notebooks/1a_filter_genomes.py

# Snakemake dry-run (notebooks only, no bio workflows)
conda run -n pyphylon-marimo snakemake -n --cores 4
```

### Docker (full pipeline execution)

```bash
# Full pipeline (notebooks + bio workflows + reports)
docker compose run pipeline snakemake --cores 4 --sdm apptainer

# Single rule
docker compose run pipeline snakemake nb_1a --cores 4

# Bio workflow only
docker compose run pipeline snakemake pangenome_mash --cores 4 --sdm apptainer

# Interactive shell inside container
docker compose run pipeline bash
```

## Pipeline DAG

```
1a → 1b → [Mash WF] → 2a → 2b → [Anno/MLST/CD-HIT WF] → 2c → 2d
                                                              ↓
3a → 3b → 4a → 5a → 5b → 5c → 5d → 5e
                                       5f ← [Infer Affinities WF]
```

- **Notebooks 1a-1b**: Genome filtering and download from BV-BRC
- **Mash WF**: Pairwise genome distance sketching
- **Notebooks 2a-2b**: Mash-based strain filtering
- **Anno/MLST/CD-HIT WF**: BAKTA annotation, MLST typing, pangenome construction
- **Notebooks 2c-2d**: Pangenome matrix generation and metadata enrichment
- **Notebooks 3a-4a**: Core/accessory/rare analysis, NMF decomposition
- **Notebooks 5a-5f**: Phylon characterization, functional enrichment, affinity inference

## Configuration

All settings live in `config.yml`:

| Key | Description |
|-----|-------------|
| `SPECIES_NAME` | Full species name (e.g., "Campylobacter jejuni") |
| `PG_NAME` | Short species code (e.g., "CJejuni") |
| `TAXON_ID` | NCBI/BV-BRC taxon ID |
| `DATA_DIR` | Path to data directory |
| `TEMP_DIR` | Path for intermediate files |
| `OUTPUT_DIR` | Path for final outputs (figures, data, reports) |
| `CONTAINER_*` | Pinned Docker image SHAs for bio tools |
| `BAKTA_THREADS` | Thread count for BAKTA annotation |

## Docker Setup

### Prerequisites

- [OrbStack](https://orbstack.dev/) (macOS) or Docker Desktop
- Docker Compose v2

### Build and Run

```bash
# Build the image (only needed when requirements.txt changes)
docker compose build

# Run the full pipeline
docker compose run pipeline snakemake --cores 4 --sdm apptainer

# Verify the environment
docker compose run pipeline python -c "import marimo, pyphylon; print('OK')"
docker compose run pipeline quarto --version
```

Code changes (notebooks, pyphylon package, Snakefile) are bind-mounted and instantly visible inside the container. Only dependency changes require a rebuild.

## Development

```bash
# Lint and format
conda run -n pyphylon-marimo ruff check notebooks/
conda run -n pyphylon-marimo ruff format notebooks/

# Run tests
pytest

# Interactive notebook development
conda run -n pyphylon-marimo marimo edit notebooks/1a_filter_genomes.py
```

## Contributing

Contributions are welcome! For major changes, please open an issue first to discuss what you would like to change. Tests use pytest.

## License

This project is licensed under the [MIT License](LICENSE).
