# Phase 1 Summary: *Campylobacter jejuni* Pangenome — Data Acquisition & Quality Filtering

**Date**: February 16–17, 2026 (PST)
**Organism**: *Campylobacter jejuni* (NCBI Taxon ID 197)
**Pipeline**: Pyphylon notebooks 1a → 1b → 2a → 2b (+ Mash workflow)
**Final dataset**: **467 genomes** (from 593 initial candidates)

---

## 1. Objective

Establish a high-quality, species-coherent genome collection for *C. jejuni* subsp. *jejuni* pangenome analysis. Phase 1 covers genome retrieval from BV-BRC, multi-stage quality filtering (assembly metrics + CheckM), genome download, metadata curation, Mash-based species verification, and hierarchical clustering.

---

## 2. Infrastructure Setup (Log: `00_setup.md`)

- Created WSL distro `pangenome` (Fedora 43) with Docker CE 29.2.1 and Miniforge conda
- Built `pyphylon:latest` Docker image (Snakemake v8.30.0 base, Debian Bookworm)
- Installed pyphylon in editable mode inside conda env `pangenome` (Python 3.11)
- Established `temp/` (intermediates) and `output/` (terminal results) directory convention
- Created `cleanup.sh` for selective phase-based cleanup of intermediate files

---

## 3. Genome Retrieval & Quality Filtering (Notebooks 1a–1b)

### 3.1 BV-BRC API Query (Log: `01_download.md`)

Queried the BV-BRC Data API for all Complete + Good *C. jejuni* genomes:
- **Endpoint**: `https://www.bv-brc.org/api/genome/`
- **Key fix**: Used `eq(taxon_lineage_ids, 197)` instead of `eq(taxon_id, 197)` to capture subspecies and strain-level taxon IDs
- **Result**: **593 genomes** retrieved with 99 metadata columns

### 3.2 Assembly Quality Filtering (Notebook 1a)

Applied three-stage quality filtering:

| Stage | Initial | Removed | Remaining |
|-------|:-------:|:-------:|:---------:|
| Pre-filtration | 593 | 0 | 593 |
| L50/N50 | 593 | 15 | 578 |
| Contig count | 578 | 0 | 578 |
| CheckM (completeness/contamination) | 578 | 85 | **493** |

- **N50 threshold**: Reference scaffold N50 retrieved via NCBI Datasets v2 REST API (replaced Selenium scraping)
- **CheckM thresholds**: Explicit Parks et al. (2015) cutoffs — completeness > 90%, contamination < 5% — used instead of kneebow auto-threshold, which failed on the uniform BV-BRC dataset
- **Notable finding**: 63 of the 85 CheckM-filtered genomes came from a single BioProject (PRJNA942088, Palestinian poultry surveillance) with elevated contamination

#### Figure: Unfiltered genome landscape (593 genomes)

Genome length vs. predicted gene count before quality filtering. Several outlier genomes visible with inflated gene counts (~2,300–2,400) relative to their genome size, indicating possible contamination.

![Unfiltered strains](../output/1a_unfiltered_strains.png)

#### Figure: Contig N50 distribution

Histogram of contig N50 values for all 593 genomes. The dashed magenta line marks the reference N50 threshold (~1.4 Mbp). Nearly all Complete genomes exceed this threshold, confirming high assembly contiguity.

![Unfiltered N50](../output/1a_unfiltered_n50.png)

#### Figure: Post-filtration genome landscape (493 genomes)

After removing 100 low-quality genomes, the genome length vs. gene count relationship is tighter and more linear, with outliers removed. Genome sizes range from ~1.5–2.0 Mbp with 1,500–2,200 predicted genes.

![Filtered strains](../output/1a_filtered_strains.png)

### 3.3 Genome Download (Notebook 1b, Log: `03_download_genome.md`)

- Downloaded 493 genome FASTA files from BV-BRC
- **Key fix**: BV-BRC FTP server now requires SSL/TLS — switched from `urllib.request.urlretrieve()` (FTP) to HTTPS Data API
- Downloads stored in `temp/1b_protected/raw/genomes/fna/` (preserved across cleanup)

### 3.4 Bug Fixes & Enhancements During Phase 1a–1b

| Issue | Fix |
|-------|-----|
| `contig_n50` assigned `contig_l50` values (line 114, `qcqa.py`) | Corrected column assignment |
| Kneebow auto-threshold too aggressive on uniform data | Added explicit `contamination_cutoff`/`completeness_cutoff` params |
| Genomes with null CheckM silently dropped | Added `checkm_missing` param (default `'keep'`) |
| Selenium scraping for reference N50 | Replaced with NCBI Datasets v2 REST API |
| FTP downloads failing | Switched to HTTPS Data API |
| numpy/numba version conflict | Pinned all dependency versions in `requirements.txt` |

---

## 4. Metadata Curation (Notebook 2a, Log: `04_clean_metadata.md`)

- Cleaned and standardized genome metadata from BV-BRC
- Fixed stale input paths (`data/interim/` → `temp/`)
- Fixed summary sync bug with `loc[]` reindexing
- Output: `temp/2a_genome_metadata.csv` and `temp/2a_genome_summary.csv` (493 genomes)

---

## 5. Mash Distance Calculation (Log: `05_run_mash.md`)

Ran all-vs-all Mash distance estimation on the 493 genomes:

- **Tool**: Mash 2.3 (installed directly in conda env; Docker bind mounts from `/mnt/f/` fail silently on WSL drvfs)
- **Sketch parameters**: Default k=21, s=1000
- **Output**: 243,049 pairwise distances (493 x 493 matrix, 28 MB TSV)
- **Sketch file**: `temp/2b_mash/combined_sketch.msh` (3.9 MB)

---

## 6. Mash Filtration & Clustering (Notebook 2b, Log: `05_run_mash.md`)

### 6.1 Species-level Mash Filtering

Used the representative strain (auto-detected medoid, genome `192222.6`) to compute pairwise distances and apply a 99th-percentile cutoff:

- **Cutoff**: 0.0257 Mash distance (~97.4% ANI)
- **Genomes removed**: 26 outliers, including:
  - 4 *C. jejuni* subsp. *doylei* genomes (subspecies-level divergence)
  - 1 likely mislabeled genome at ~90.7% ANI
  - Other borderline/divergent strains
- **Remaining**: **467 genomes**

#### Figure: Pre-filtration Mash distance distribution

Histogram of all 243,049 pairwise Mash distances. The bulk of distances fall between 0.005–0.030 (within-species), with a long tail extending to ~0.12 representing subspecies and mislabeled genomes.

![Mash distances — all](../output/2b_mash_dist_all.png)

#### Figure: Post-filtration Mash distance distribution (467 genomes)

After removing outliers, distances are concentrated between 0.000–0.030 with no long tail. The multimodal structure (peaks near 0.000, 0.013, and 0.020) reflects subpopulation structure within *C. jejuni* subsp. *jejuni*.

![Mash distances — filtered](../output/2b_mash_dist_filtered.png)

#### Figure: Mash distance heatmap

Full 493 x 493 Mash distance matrix. Bright (yellow) cells indicate high distance. A clear block-diagonal structure is visible, with a handful of highly divergent genomes forming bright rows/columns — these are the outliers removed by the 99th-percentile filter.

![Mash heatmap](../output/2b_mash_heatmap.png)

### 6.2 Hierarchical Clustering & Sensitivity Analysis

Applied hierarchical clustering on the Pearson correlation distance matrix:

- **Sensitivity analysis**: Swept distance thresholds from 0 to 1 and tracked cluster count
- **Elbow point**: 92 clusters at threshold 0.30 (cluster count decelerates beyond this point)

#### Figure: Sensitivity analysis — cluster count vs. threshold

The curve shows exponential decay of cluster count as the distance threshold increases. The dashed line marks the elbow at 92 clusters (threshold = 0.30), where the rate of cluster merging decelerates.

![Sensitivity analysis](../output/2b_sensitivity_analysis.png)

#### Figure: Hierarchical clustermap (initial, 20 clusters)

Seaborn clustermap of the 467-genome correlation distance matrix with hierarchical dendrogram and 20-cluster color annotation. Clear block-diagonal structure with several well-defined clades visible.

![Clustermap — initial](../output/2b_clustermap_initial.png)

#### Figure: Cluster size distribution (initial)

Right-skewed distribution with most clusters containing 1–10 genomes and a few large clusters (30–55 genomes). This reflects the known clonal complex structure of *C. jejuni*.

![Cluster sizes — initial](../output/2b_cluster_sizes_initial.png)

#### Figure: Hierarchical clustermap (final, 20 clusters)

After stability filtering (small_clst_limit = 0, no additional removal), the final clustermap is identical to the initial. All 467 genomes retained.

![Clustermap — final](../output/2b_clustermap_final.png)

#### Figure: Cluster size distribution (final)

Final cluster size distribution. Slightly more concentrated than initial, with the largest clusters containing ~55 genomes.

![Cluster sizes — final](../output/2b_cluster_sizes_final.png)

---

## 7. Literature Collection (Logs: `03_literature_search.md`, `06_paper_download.md`)

In parallel with the computational pipeline:

- Curated **52 references** on *C. jejuni* pangenomics, AMR, virulence, MLST, and population structure
- Covered 14 BioProject accessions used in existing studies
- Automated PDF download: retrieved **46/52 papers** (84 MB) via Europe PMC, Unpaywall, and arXiv APIs
- References stored in `papers/papers.md` with `[PDF]` links

---

## 8. Pipeline Preparation for Phase 2 (Log: `07_build_cds.md`)

Updated notebook `2c_build_cds_pangenome.ipynb` for the next phase:

- Fixed paths to use `temp/` + `output/` convention
- Added `SNAKEMAKE_DATA_DIR` to config.yml to separate Snakemake outputs
- Rewrote genome filtering with O(n) set-based lookup (from O(n^2) nested loop)
- Pipeline is ready for BAKTA annotation and CD-HIT pangenome construction

---

## 9. Summary of Genome Counts

| Pipeline Stage | Notebook | Genomes | Removed |
|----------------|:--------:|:-------:|:-------:|
| BV-BRC query (Complete + Good) | 1a | 593 | — |
| L50/N50 filtering | 1a | 578 | 15 |
| CheckM completeness/contamination | 1a | 493 | 85 |
| Genome download | 1b | 493 | 0 |
| Metadata curation | 2a | 493 | 0 |
| Mash species filtering (99th %ile) | 2b | 467 | 26 |
| Cluster stability filtering | 2b | 467 | 0 |
| **Final dataset** | — | **467** | **126 total** |

---

## 10. Output Files

### Terminal Results (`output/`)

| File | Description |
|------|-------------|
| `1a_df_filtration.csv` | Filtration report (genome counts per stage) |
| `1a_unfiltered_strains.png` | Genome length vs. gene count (pre-filter) |
| `1a_unfiltered_n50.png` | Contig N50 distribution with reference threshold |
| `1a_filtered_strains.png` | Genome length vs. gene count (post-filter) |
| `2b_mash_square.csv` | 467 x 467 raw Mash distance matrix |
| `2b_mash_corr_dist.csv` | 467 x 467 correlation distance matrix |
| `2b_genome_summary.csv` | Final 467-genome summary |
| `2b_mash_heatmap.png` | Full Mash distance heatmap |
| `2b_mash_dist_all.png` | Distance histogram (pre-filtration) |
| `2b_mash_dist_filtered.png` | Distance histogram (post-filtration) |
| `2b_sensitivity_analysis.png` | Cluster count vs. threshold elbow plot |
| `2b_clustermap_initial.png` | Hierarchical clustermap (initial) |
| `2b_clustermap_final.png` | Hierarchical clustermap (final) |
| `2b_cluster_sizes_initial.png` | Cluster size distribution (initial) |
| `2b_cluster_sizes_final.png` | Cluster size distribution (final) |

### Intermediate Files (`temp/`)

| File | Description |
|------|-------------|
| `1a_genome_summary.csv`, `1a_genome_metadata.csv` | Post-QC genome data |
| `1b_genome_summary.csv`, `1b_genome_metadata.csv` | Post-download genome data |
| `1b_protected/raw/genomes/fna/*.fna` | 493 downloaded genome FASTA files |
| `2a_genome_summary.csv`, `2a_genome_metadata.csv` | Post-metadata-curation data |
| `2b_genome_metadata.csv` | Final filtered metadata (consumed by Snakemake) |
| `2b_mash/combined_sketch.msh` | Mash sketch file (3.9 MB) |

---

## 11. Next Steps (Phase 2)

1. **BAKTA annotation**: Run genome annotation workflow via Snakemake on 467 genomes
2. **MLST typing**: Determine multi-locus sequence types for all strains
3. **CD-HIT pangenome construction**: Cluster coding sequences to build the strain-by-gene presence/absence matrix (P matrix)
4. **Notebook 2c–2d**: Generate pangenome matrices and enrich with metadata
