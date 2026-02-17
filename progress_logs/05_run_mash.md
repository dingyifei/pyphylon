# 05 — Run Mash Workflow

**Date**: 2026-02-17 ~4:19 AM PST

## Summary

Ran Mash sketch + distance on 493 *C. jejuni* genomes to produce pairwise distance matrix for notebook 2b.

## Steps

1. **Checked Docker** — Docker CE already running in WSL `pangenome` distro.
2. **Docker mount issue** — Native Docker CE in WSL cannot mount `/mnt/f/` (Windows drvfs) paths into containers. Bind mounts create empty directories.
3. **Installed mash directly** — `conda install -y -c bioconda mash` in `pangenome` conda env (mash 2.3, plus capnproto 1.0.2 and gsl 2.8 dependencies).
4. **Created output directory** — `examples/data/processed/mash/`
5. **Ran `mash sketch`** — Sketched all 493 genomes from `temp/1b_protected/raw/genomes/fna/*.fna` into `combined_sketch.msh` (3.9 MB).
6. **Ran `mash dist`** — All-vs-all distance calculation → `mash_distances.txt` (28 MB, 243,049 rows = 493×493).
7. **Verified** — Row count matches expected, file format has 5 TSV columns (genome1, genome2, mash_distance, p_value, matching_hashes). Genome paths use full `temp/1b_protected/...` prefix but notebook 2b extracts just the genome ID via `x.split('/')[-1].split('.fna')[0]`.

## Output Files

| File | Size | Description |
|------|------|-------------|
| `examples/data/processed/mash/combined_sketch.msh` | 3.9 MB | Mash sketch of 493 genomes |
| `examples/data/processed/mash/mash_distances.txt` | 28 MB | Pairwise distance matrix (243,049 rows) |

## Lessons Learned

- **Docker bind mounts in WSL**: Native Docker CE running inside a WSL distro cannot mount paths from `/mnt/f/` (Windows drives via drvfs). The mount succeeds without error but the container sees an empty directory. Workaround: install tools directly in WSL and run without Docker, or copy data to native Linux filesystem first.
- **Mash via conda**: `mash` is available in bioconda channel, installs cleanly in seconds. No need for Docker/Singularity for this simple two-command workflow.

## Access

```bash
# View distance matrix
head examples/data/processed/mash/mash_distances.txt

# Mash is now available in WSL
wsl -d pangenome bash -c "source ~/miniforge3/etc/profile.d/conda.sh && conda activate pangenome && mash --version"
```

## Update: Fix Hardcoded Representative Strain in Notebook 2b

**Date**: 2026-02-17 ~5:03 AM PST

Cell 18 had a hardcoded representative strain ID `'1314.132'` from the original *S. pyogenes* example, causing a `KeyError` when running with *C. jejuni* data.

**Fix** (cell 18): Replaced hardcoded ID with auto-detection logic:
1. Looks up reference genomes from the `reference_genome` column in `scrubbed_species_metadata`.
2. Filters to only those present in the mash distance matrix (`df_mash_square.index`).
3. Falls back to the **medoid** (genome with smallest mean mash distance) if no reference is annotated.

This makes the notebook organism-agnostic — for *C. jejuni* it resolves to `['192222.6']`.

## What's Next

Notebook 2b (`2b_mash_filtration_and_clustering.ipynb`) can now be run — it reads `data/processed/mash/mash_distances.txt` and produces Mash-filtered metadata/summary CSVs in `output/`.
