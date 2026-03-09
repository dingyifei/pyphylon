# 08 — BAKTA Annotation, MLST, and CD-HIT Workflow

**Date**: 2026-02-17 06:40 PST

## Overview

Run the BAKTA annotation, MLST typing, and CD-HIT pangenome clustering pipeline on 467 filtered *C. jejuni* genomes via Docker + Snakemake on native Linux FS (workaround for drvfs bind mount issue).

## Key Steps

1. **Fixed `build_cds_pangenome()` in `pyphylon/pangenome.py`** (lines 106-109)
   - Added `.cdhit` suffix to `output_nr_faa_copy` (was empty string)
   - Fixed `.cdhit.clstr` suffix on `output_nr_clstr` (was `.clstr`, with typo "clusteeipxzr")
   - Uncommented `cluster_with_cdhit()` call and `os.remove()` cleanup
   - Now matches the working `build_noncoding_pangenome()` pattern (lines 368-371)

2. **Created native Linux FS directory structure**
   ```
   /home/yifei/pyphylon_run/
   ├── examples/
   │   ├── config.yml
   │   ├── data/raw/genomes/fna/  (493 .fna files)
   │   ├── data/processed/        (empty, output target)
   │   └── temp/2b_genome_metadata.csv
   └── workflow/anno_mlst_cdhit/Snakefile
   ```

3. **Copied data from Windows FS to native Linux FS**
   - 493 genome files (~510 MB) from `temp/1b_protected/raw/genomes/fna/`
   - `config.yml`, `2b_genome_metadata.csv`, `Snakefile`
   - Reason: Docker bind mounts from `/mnt/f/` (drvfs) fail silently

4. **Started Docker daemon** in WSL (`sudo dockerd &`)

5. **Launched Snakemake in detached Docker container**
   ```bash
   docker run --privileged -d --name pyphylon_anno \
     -v /home/yifei/pyphylon_run/examples:/examples \
     -v /home/yifei/pyphylon_run/workflow:/workflow \
     pyphylon:latest \
     bash -c "cd /workflow/anno_mlst_cdhit && snakemake -d /examples/data --use-singularity -c 32 2>&1 | tee /examples/snakemake.log"
   ```
   - Container ID: `6faa21536ab6`
   - Container name: `pyphylon_anno`
   - Using `-c 32` (all cores) for maximum parallelism

## Workflow Details

The Snakefile runs these rules in order:
1. **get_db** — Downloads BAKTA db-light (~1.2 GB from Zenodo, one-time)
2. **bakta_annotation** — Annotates all 467 genomes (8 threads/job → 4 parallel with -c 32)
3. **mlst** — Types all 467 genomes (4 threads/job → 8 parallel)
4. **concat_faas** + **cdhit** — Concatenates .faa files, clusters at 80% identity
5. **mlst_report** — Concatenates all MLST results

**Estimated runtime**: ~10 hours with `-c 32` (4 parallel BAKTA jobs × ~5 min each for 467 genomes)

## Access Instructions

### Monitor progress
```bash
# From Git Bash:
MSYS_NO_PATHCONV=1 wsl -d pangenome bash -c 'docker logs pyphylon_anno --tail 50'

# Or check the log file:
MSYS_NO_PATHCONV=1 wsl -d pangenome bash -c 'tail -50 /home/yifei/pyphylon_run/examples/snakemake.log'
```

### Check container status
```bash
MSYS_NO_PATHCONV=1 wsl -d pangenome bash -c 'docker ps --filter name=pyphylon_anno'
```

### If interrupted — resume
```bash
# Re-enter container:
MSYS_NO_PATHCONV=1 wsl -d pangenome bash -c 'docker start pyphylon_anno && docker exec -it pyphylon_anno bash'

# Inside container:
cd /workflow/anno_mlst_cdhit
snakemake -d /examples/data --use-singularity -c 32
```

### After completion — copy results back
```bash
MSYS_NO_PATHCONV=1 wsl -d pangenome bash -c '
  PROJECT=/mnt/f/lab_projects/pangenomics/pyphylon
  RUN=/home/yifei/pyphylon_run
  mkdir -p $PROJECT/examples/data/processed
  cp -r $RUN/examples/data/processed/bakta $PROJECT/examples/data/processed/
  cp -r $RUN/examples/data/processed/mlst $PROJECT/examples/data/processed/
  cp $RUN/examples/data/processed/mlst_report.txt $PROJECT/examples/data/processed/
  cp -r $RUN/examples/data/processed/cd-hit-results $PROJECT/examples/data/processed/
'
```

### Verify results
```bash
# BAKTA: should have 467 subdirectories, each with .faa/.gbff/.tsv
ls examples/data/processed/bakta/ | wc -l  # expect 467

# MLST: should have 467 .tsv files + mlst_report.txt
ls examples/data/processed/mlst/ | wc -l   # expect 467

# CD-HIT: should have CJejuni file + .clstr
ls examples/data/processed/cd-hit-results/  # expect CJejuni, CJejuni.clstr
```

6. **Updated test expectations** in `pyphylon/test/test_pangenome.py`
   - With CD-HIT now running, allele count changed from 5279 → 4218 (clustering reduces unique sequences)
   - All 4 pangenome tests pass
   - Installed `cd-hit` 4.8.1 in pangenome conda env via bioconda

## Status

- [x] pangenome.py fix applied
- [x] Test expectations updated (4/4 pass)
- [x] cd-hit installed in pangenome conda env
- [x] Data copied to native Linux FS
- [x] Docker container launched
- [x] Snakemake workflow complete (939/939 steps, ~2h 3min runtime)
- [x] Results copied back to Windows FS
- [x] Verification complete (467 BAKTA dirs, 467 MLST files, CD-HIT CJejuni + .clstr)
