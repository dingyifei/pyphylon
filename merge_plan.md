# Merge Plan: Quarto Reports → Marimo Notebooks

> **Goal:** Merge narrative text from Quarto `.qmd` reports and original Jupyter notebook comments into marimo notebooks as `mo.md()` cells. Add light interactive features. Then remove the Quarto layer entirely.
>
> **How to resume:** _"Read merge_plan.md and execute the next unchecked stage."_

---

## References

| Source | Path | Purpose |
|--------|------|---------|
| Quarto reports | `reports/*.qmd` | Polished narrative text, section structure |
| Original Jupyter notebooks | `/mnt/f/lab_projects/pangenomics/ref_repo/pyphylon/examples/` | Markdown cells, inline comments, original flow |
| Current marimo notebooks | `notebooks/*.py` | Target files to modify |
| Snakefile | `Snakefile` (lines 272–462) | Report rules to remove |

---

## Conventions

### Markdown cells (narrative-only)

Return `mo.md()` as the last expression — NOT wrapped in `mo.output.replace()`. This makes marimo auto-detect them as markdown cells (special editor) and populates the sidebar outline.

```python
@app.cell
def _():
    mo.md(
        """
        # 1a: Filter Genomes for Download

        Query genomes from BV-BRC, filter by species, then quality-filter
        using N50, contig count, and CheckM metrics.
        """
    )
```

- **H1** (`#`) for notebook title — one per notebook
- **H2** (`##`) for sections (Setup, Plot, Filter, Save, etc.)
- **H3** (`###`) for subsections

### Computation cells (dynamic output)

Keep existing `mo.output.replace(mo.md(...))` or `mo.output.replace(mo.vstack([...]))` pattern for cells that show runtime values.

### Interactive features (low-effort additions)

Apply these where they fit naturally — don't force them everywhere:

| Feature | When to use | Example |
|---------|------------|---------|
| `mo.ui.table(df)` | Any DataFrame display (sorting, filtering built-in) | Filtration reports, enrichment results, summary stats |
| `mo.ui.plotly(fig)` | Plotly figures (hover, zoom, pan) | Circular genome plots (already Plotly in `plotting.py`) |
| `mo.ui.tabs({...})` | Multi-phylon outputs | Wordclouds, circular plots, per-phylon enrichments |
| `mo.accordion({...})` | Optional detail sections | Summary stats, debug info |
| `mo.lazy(...)` | Expensive renders inside tabs/accordion | Wrap heavy plots so they only compute when opened |

### Commit convention

One commit per stage: `Merge report text into nb_<ID>` (e.g., `Merge report text into nb_1a`)

---

## Checklist

### Stage 1: Notebook 1a — Filter Genomes
- [x] **Marimo:** `notebooks/1a_filter_genomes.py`
- [x] **Quarto:** `reports/1a_genome_filtering.qmd`
- [x] **Ref Jupyter:** `ref_repo/pyphylon/examples/1a_filter_genomes_for_download.ipynb`
- [x] Add H1 title + overview cell (before config cell)
- [x] Add `## Setup` header cell (before config cell)
- [x] Add `## Plot Unfiltered Dataset` header cell (before unfiltered plot cell)
- [x] Add `## Quality Filtering` header cell (before filter cell)
- [x] Add `## Filtered Distributions` header cell (before filtered plot cell)
- [x] Add `## Save Filtered Genomes` header cell (before save cell)
- [x] Ensure `mo.ui.table(df_filtration)` is kept (already there)
- [x] Lint: `ruff check notebooks/1a_filter_genomes.py`
- [x] Commit

### Stage 2: Notebook 1b — Download Genomes
- [x] **Marimo:** `notebooks/1b_download_genomes.py`
- [x] **Quarto:** `reports/1b_genome_download.qmd`
- [x] **Ref Jupyter:** `ref_repo/pyphylon/examples/1b_download_genomes_bvbrc.ipynb`
- [x] Add H1 title + overview cell
- [x] Add section headers matching ref notebook flow
- [x] Lint + commit

### Stage 3: Notebook 2a — Clean Metadata
- [x] **Marimo:** `notebooks/2a_clean_metadata.py`
- [x] **Quarto:** `reports/2a_clean_metadata.qmd`
- [x] **Ref Jupyter:** `ref_repo/pyphylon/examples/2a_clean_metadata.ipynb`
- [x] Add H1 title + overview (de-duplication by biosample_accession)
- [x] Add section headers matching ref notebook flow
- [x] Lint + commit

### Stage 4: Notebook 2b — Mash Filtration
- [ ] **Marimo:** `notebooks/2b_mash_filtration.py`
- [ ] **Quarto:** `reports/2b_mash_filtration.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/2b_mash_filtration_and_clustering.ipynb`
- [ ] Add H1 title + overview (ANI-based filtering, hierarchical clustering)
- [ ] Add `## Raw Mash Distances` before heatmap cell
- [ ] Add `## Distance Distribution` before histogram cells
- [ ] Add `## Sensitivity Analysis` before sensitivity cell
- [ ] Add `## Initial Clustering` before initial clustermap cell
- [ ] Add `## Final Clustering` before final clustermap cell
- [ ] Add `## Summary Statistics` before summary cell
- [ ] Consider `mo.ui.table()` for cluster size summary
- [ ] Lint + commit

### Stage 5: Notebook 2c — Build Pangenome
- [ ] **Marimo:** `notebooks/2c_build_pangenome.py`
- [ ] **Quarto:** `reports/2c_build_pangenome.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/2c_build_cds_pangenome.ipynb`
- [ ] Add H1 title + overview (CD-HIT clustering → gene/allele matrices)
- [ ] Add `## Genes per Genome Distribution` before plot cell
- [ ] Lint + commit

### Stage 6: Notebook 2d — Enrich Metadata
- [ ] **Marimo:** `notebooks/2d_enrich_metadata.py`
- [ ] **Quarto:** `reports/2d_enrich_metadata.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/2d_enrich_metadata.ipynb`
- [ ] Add H1 title + overview (MLST join)
- [ ] Add `## MLST Sequence Type Distribution` before MLST cell
- [ ] Consider `mo.ui.table()` for MLST counts
- [ ] Lint + commit

### Stage 7: Notebook 3a — Extract CAR
- [ ] **Marimo:** `notebooks/3a_extract_CAR.py`
- [ ] **Quarto:** `reports/3a_extract_CAR.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/3a_extract_CAR_genomes.ipynb`
- [ ] Add H1 title + overview (Core/Accessory/Rare, power-law model)
- [ ] Add `## Gene Frequency Distribution` before frequency plot
- [ ] Add `## Pangenome Segment Fitting` before segments plot
- [ ] Consider `mo.ui.table()` for CAR summary
- [ ] Lint + commit

### Stage 8: Notebook 3b — Heaps Plot
- [ ] **Marimo:** `notebooks/3b_heaps_plot.py`
- [ ] **Quarto:** `reports/3b_heaps_plot.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/3b_heaps_plot.ipynb`
- [ ] Add H1 title + overview with Heaps' formula: $G(n) = \kappa \cdot n^\gamma$
- [ ] Add note: gamma < 1 → open pangenome
- [ ] Add `## Heaps' Law Stacked Plot` before plot cell
- [ ] Lint + commit

### Stage 9: Notebook 4a — NMF Decomposition
- [ ] **Marimo:** `notebooks/4a_nmf_decomposition.py`
- [ ] **Quarto:** `reports/4a_nmf_decomposition.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/4a_NMF_decomposition.ipynb`
- [ ] Add H1 title + overview (NMF, MCA rank selection, consensus clustering)
- [ ] Add `## MCA Cumulative Variance` before MCA plot
- [ ] Add `## Consensus Matrix` before consensus clustermap
- [ ] Add `## Consensus Matrix (Filtered)` + explanation (excludes 50–100%, 75–100% submatrices)
- [ ] Add `## NMF Output Shapes` before summary
- [ ] Consider `mo.ui.table()` for NMF summary
- [ ] Lint + commit

### Stage 10: Notebook 5a — Phylon Characterization
- [ ] **Marimo:** `notebooks/5a_phylon_characterization.py`
- [ ] **Quarto:** `reports/5a_phylon_characterization.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/5a_phylon_characterization.ipynb`
- [ ] Add H1 title + overview (L binarized via 3-means, A binarized via threshold)
- [ ] Add `## Gene Frequency vs Phylon Count` + regression explanation
- [ ] Add `## Sorted L_binarized` + sorting explanation (zero/single/poly-phylon)
- [ ] Add `## Sorted A_binarized` + sorting explanation
- [ ] Wrap circular genome plots with `mo.ui.plotly()` if applicable
- [ ] Lint + commit

### Stage 11: Notebook 5b — Gene Differentiation
- [ ] **Marimo:** `notebooks/5b_gene_diff.py`
- [ ] **Quarto:** `reports/5b_gene_diff.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/5b_gene_diff_top_down.ipynb`
- [ ] Add H1 title + overview (Ward's method, shared gene content)
- [ ] Add `## Ward Clustering Heatmap` before clustermap
- [ ] Add `## Phylon Dendrogram` + annotation explanation (exclusive genes, split genes)
- [ ] Lint + commit

### Stage 12: Notebook 5c — Functional Enrichments
- [ ] **Marimo:** `notebooks/5c_functional_enrichments.py`
- [ ] **Quarto:** `reports/5c_functional_enrichments.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/5c_functional_enrichments.ipynb`
- [ ] Add H1 title + overview (GO enrichment, hypergeometric test, p < 0.05)
- [ ] Add `## Top Enrichments` — use `mo.ui.table()` for enrichment results
- [ ] Add `## Enrichment Heatmap` + explanation (p-values, darker = stronger)
- [ ] Add `## Phylon Wordclouds` — organize with `mo.ui.tabs()` per phylon
- [ ] Lint + commit

### Stage 13: Notebook 5d — Gene Alignment
- [ ] **Marimo:** `notebooks/5d_gene_alignment.py`
- [ ] **Quarto:** `reports/5d_gene_alignment.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/5d_gene_alignment_visualization.ipynb`
- [ ] Add H1 title + overview (anchor gene alignment, gene order standardization)
- [ ] Add `## Genetic Variation` before variation cell
- [ ] Add `## Gene Length Distribution` before histogram
- [ ] Add `## Circular Genome Plots` + explanation — organize with `mo.ui.tabs()` per phylon
- [ ] Wrap circular plots with `mo.ui.plotly()` for interactive hover
- [ ] Add `## Unique Genes by Phylon` — use `mo.ui.table()`
- [ ] Lint + commit

### Stage 14: Notebook 5e — BLAST Enrichment
- [ ] **Marimo:** `notebooks/5e_blast_enrichment.py`
- [ ] **Quarto:** `reports/5e_blast_enrichment.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/5e_pangenome_local_blastdb_enrichment.ipynb`
- [ ] Add H1 title + overview (VFDB comparison)
- [ ] Add `## VFDB Enrichment Results` — use `mo.ui.table()` for hits
- [ ] Add `## Custom Queries` with blastp command example
- [ ] Lint + commit

### Stage 15: Notebook 5f — Infer Affinities
- [ ] **Marimo:** `notebooks/5f_infer_affinities.py`
- [ ] **Quarto:** `reports/5f_infer_affinities.qmd`
- [ ] **Ref Jupyter:** `ref_repo/pyphylon/examples/5f_infer_affinities.ipynb`
- [ ] Add H1 title + overview (NNLS affinity estimation, 5-NN MASH cluster classification)
- [ ] Add `## MASH Cluster Assignments` — use `mo.ui.table()`
- [ ] Add `## Affinity Heatmap`
- [ ] Lint + commit

### Stage 16: Snakefile Cleanup
- [ ] Remove all 15 `report_*` rules (lines 272–462)
- [ ] Remove `REP = config["REPORTS_DIR"]` (line 14)
- [ ] Update `rule all:` — replace report PDF targets with notebook data outputs
- [ ] Update header comment (remove Quarto references)
- [ ] Dry-run: `snakemake -n --cores 1`
- [ ] Commit: `Remove Quarto report rules from Snakefile`

### Stage 17: Delete Quarto Files
- [ ] Delete `reports/*.qmd` (15 files)
- [ ] Delete `reports/*.quarto_ipynb` (3 files)
- [ ] Delete `reports/_quarto.yml`, `reports/.gitignore`
- [ ] Delete `_quarto.yml` (root)
- [ ] Remove `reports/` directory
- [ ] Remove `REPORTS_DIR` from `config.yml` (line 10)
- [ ] Commit: `Remove Quarto report files and config`

### Stage 18: Documentation Update
- [ ] Update `CLAUDE.md` — remove Quarto references, `REPORTS_DIR`, `report_*` examples
- [ ] Verify: `grep -r "quarto\|REPORTS_DIR\|report_" Snakefile config.yml CLAUDE.md`
- [ ] Delete this file (`merge_plan.md`)
- [ ] Commit: `Update docs, remove merge plan`

---

## Progress Log

| Date | Stage | Commit | Notes |
|------|-------|--------|-------|
| 2026-03-03 | 1 — nb_1a | `Merge report text into nb_1a` | 6 markdown cells added (H1, Setup, Plot Unfiltered, Quality Filtering, Filtered Distributions, Save) |
| 2026-03-03 | 2 — nb_1b | `Merge report text into nb_1b` | 5 markdown cells added (H1, Setup, Download, Update Genome Files, Save); removed unused DEBUG variable |
| 2026-03-03 | 3 — nb_2a | `Merge report text into nb_2a` | 4 markdown cells added (H1 + overview, Setup, De-duplicate Entries, Save Files) |
