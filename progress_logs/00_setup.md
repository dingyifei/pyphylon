# 00 - WSL Environment Setup

**Date**: 2026-02-16 ~19:43–20:15 PST

## Summary

Set up a project-local WSL environment named `pangenome` based on Fedora 43, with user account, Docker, conda environment, and pyphylon installed.

## Key Steps

1. Listed existing WSL distributions — identified `FedoraLinux-43` as the latest Fedora distro
2. Created `.wsl\` directory under the project root
3. Exported `FedoraLinux-43` to `.wsl\fedora43.tar`, imported as `pangenome`, removed tar
4. Created user `yifei` with passwordless sudo, set as default user via `/etc/wsl.conf`
5. Installed Docker CE 29.2.1 (`dnf install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin`), added `yifei` to docker group
6. Installed Miniforge (conda) at `/home/yifei/miniforge3`, initialized for bash
7. Created conda environment `pangenome` (Python 3.11) with all dependencies from `requirements.txt` + `conda/environment.yml`
8. Installed pyphylon in editable mode (`pip install -e .`) inside the conda env
9. Updated Dockerfile base image from `snakemake/snakemake:v7.25.0` (Debian Buster, EOL) to `snakemake/snakemake:v8.30.0` (Debian Bookworm, micromamba-based). Replaced `apt install` with `micromamba install` for `less`.
10. Built Docker image `pyphylon:latest` (3.96GB, based on snakemake/snakemake:v8.30.0)
11. Verified: all Python imports pass, Docker image present, snakemake 8.30.0 runs, `--use-singularity` flag still supported

## Access

### WSL Distro
```bash
# Launch the environment (defaults to user yifei)
wsl -d pangenome

# From Git Bash (avoid path mangling)
MSYS_NO_PATHCONV=1 wsl -d pangenome -- <command>
```

### Conda Environment
```bash
# Inside the WSL
conda activate pangenome
python -c "import pyphylon; print('OK')"
```

### Docker
```bash
# Start Docker daemon (if not running)
sudo dockerd &
# or
sudo systemctl start docker

# Run pyphylon container with workflows
docker run --privileged -it \
  -v $(pwd)/examples:/examples \
  -v $(pwd)/workflow:/workflow \
  pyphylon
```

## Environment Details

| Component | Details |
|-----------|---------|
| WSL distro | `pangenome` |
| Base OS | Fedora 43 |
| WSL version | 2 |
| Disk location | `F:\lab_projects\pangenomics\pyphylon\.wsl\pangenome\` |
| Default user | `yifei` (passwordless sudo) |
| Docker | CE 29.2.1 |
| Conda | Miniforge at `/home/yifei/miniforge3` |
| Conda env | `pangenome` (Python 3.11.14) |
| Docker image | `pyphylon:latest` (3.96GB, snakemake v8.30.0) |

## Notes

- Dockerfile base image upgraded from v7.25.0 to v8.30.0 — this change is local and uncommitted
- v8.30.0 uses micromamba instead of apt/dpkg; packages installed via `micromamba install`
- `--use-singularity` still works in v8 (aliased to `--use-apptainer`)
- All v7.x images use Debian Buster (EOL); only v8+ moved to Bookworm
- Docker daemon needs to be started manually in WSL (no systemd auto-start by default)
