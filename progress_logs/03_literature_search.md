# Literature Search — C. jejuni Pangenome & Genomics

**Date**: 2026-02-17 ~2:57 AM – 3:16 AM PST

## Summary

Created `papers.md` at project root with curated C. jejuni pangenome and genomics literature.

## Key Steps

1. Extracted 14 unique BioProject accessions from `output/1a_genome_metadata.csv`
2. Ran Scholar Gateway semantic searches (2 calls) for C. jejuni pangenome/population structure and AMR/virulence literature
3. Ran WebSearch queries (10+ calls) covering:
   - Pangenome studies
   - Core/accessory genome & population structure
   - Antimicrobial resistance (AMR)
   - Virulence & pathogenicity
   - MLST & phylogenomics
   - NCTC 11168 reference genome
   - NMF/pangenome matrix decomposition methodology
4. Fetched all 14 NCBI BioProject pages via WebFetch to identify associated publications and submitter organizations
5. Spot-checked 5 PubMed links and 3 DOI redirects for accuracy
6. Fixed 4 author/PMID errors found during verification (Zhong not Gupta, Chauhan not Monk, Peters not Konkel, Skarp not de Haan + wrong PMID)
7. Confirmed 0 duplicate DOIs across all 46 entries
8. Hard-coded 10 methodology references (Roary, Panaroo, CD-HIT, Mash, BAKTA, CheckM, UMAP, HDBSCAN, NMF decomposition, NMF review)

## Output

- **File**: `papers.md` (project root, 219 lines)
- **Total papers**: 52 entries (46 with DOIs, 6 BioProject entries noting "no publication")
- **Sections**: 8 (Pangenome, Core/Accessory, AMR, Virulence, MLST, Reference Genome, Methodology, BioProject-Associated)
- **BioProjects covered**: All 14 from metadata
  - 3 with publications found (PRJNA168989: 2 papers, PRJNA1087001: 1 paper, PRJNA492384: 1 related paper)
  - 10 marked as direct NCBI submissions (no dedicated publication)
  - 1 flagged as legacy/metadata error (PRJNA175 → Prevotella intermedia)

## Notes

- PubMed MCP server was persistently rate-limited after initial burst; all searches completed via WebSearch and Scholar Gateway instead
- Scholar Gateway results were mostly Wiley-published papers, so coverage skewed toward Environmental Microbiology and Molecular Ecology journals
- BioProjects are mostly recent (2024–2025) direct submissions without publications yet
