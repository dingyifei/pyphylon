import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import pandas as pd
    import yaml

    from pyphylon.blast_utils import (
        blast_localdb_enrichment,
        extract_reference_dna_sequences,
        extract_reference_sequences,
        make_blast_db,
        process_blast_results,
    )


@app.cell
def _():
    mo.md("""
    # 5e: BLAST Enrichment Analysis

    Compare pangenome representative sequences against external databases
    to identify functionally relevant genes. This notebook:

    1. **Extracts representative alleles** — one protein and one DNA sequence
       per CD-HIT gene cluster
    2. **Searches VFDB** — BLASTs representative proteins against the
       [Virulence Factor Database](https://www.mgc.ac.cn/VFs/download.htm)
       (e-value < 1e-5, identity > 80%)
    3. **Builds a pangenome BLAST DB** — enables custom blastp queries
       against the full pangenome
    """)
    return


@app.cell
def _():
    mo.md("""
    ## Setup
    """)
    return


@app.cell
def _():
    """Parse config and set up directories."""
    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    DATA = CONFIG["DATA_DIR"]
    TEMP = CONFIG["TEMP_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]
    SPECIES = CONFIG["PG_NAME"]

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    return DATA, OUT, SPECIES


@app.cell
def _():
    mo.md("""
    ## Extract Representative Sequences

    Extract one representative protein and one DNA sequence per CD-HIT gene
    cluster. These are the longest member of each cluster and serve as the
    query set for downstream BLAST searches.
    """)
    return


@app.cell
def _(DATA, SPECIES):
    """Extract representative protein and DNA sequences from CD-HIT results."""
    cd_hit_dir = os.path.join(DATA, "processed/cd-hit-results")
    prot_seqs = os.path.join(cd_hit_dir, f"{SPECIES}_representative_sequences")
    dna_seqs = os.path.join(cd_hit_dir, f"{SPECIES}_representative_DNA_sequences")

    if not os.path.exists(prot_seqs):
        extract_reference_sequences(cd_hit_dir, SPECIES, prot_seqs)
        mo.output.replace(mo.md(f"Extracted representative protein sequences to `{prot_seqs}`"))
    else:
        mo.output.replace(mo.md(f"Representative protein sequences already exist at `{prot_seqs}`"))

    if not os.path.exists(dna_seqs):
        extract_reference_dna_sequences(DATA, SPECIES, dna_seqs)
    return dna_seqs, prot_seqs


@app.cell
def _():
    mo.md("""
    ## VFDB Enrichment Results

    Pangenome representative protein sequences are searched against
    [VFDB](https://www.mgc.ac.cn/VFs/download.htm) (Virulence Factor
    Database) using blastp with e-value < 1e-5 and identity > 80%.

    > **Prerequisite:** Download `VFDB_setA_pro.fas` from
    > <https://www.mgc.ac.cn/VFs/download.htm> and place it in
    > `data/external/VFDB/`.
    """)
    return


@app.cell
def _(DATA, prot_seqs):
    """BLAST pangenome representative sequences against VFDB (optional).

    Requires manual download of VFDB_setA_pro.fas from
    https://www.mgc.ac.cn/VFs/download.htm into data/external/VFDB/.
    """
    vfdb_fasta = os.path.join(DATA, "external/VFDB/VFDB_setA_pro.fas")
    vfdb_db = os.path.join(DATA, "external/VFDB/VFDB")
    vfdb_results_path = os.path.join(DATA, "external/VFDB/results.txt")

    vfdb_results = None
    if os.path.exists(vfdb_fasta):
        if not os.path.exists(vfdb_db + ".phr"):
            make_blast_db(vfdb_fasta, vfdb_db)
        if not os.path.exists(vfdb_results_path):
            blast_localdb_enrichment(vfdb_db, prot_seqs, vfdb_results_path, e_val=1e-5)
        vfdb_results = process_blast_results(vfdb_results_path, e_val=1e-5, percent_identity=80)
        mo.output.replace(
            mo.vstack([
                mo.md(
                    f"VFDB BLAST: **{len(vfdb_results)}** hits "
                    f"({vfdb_results['query'].nunique()} unique query genes, "
                    f"{vfdb_results['target'].nunique()} unique VFDB targets)"
                ),
                mo.ui.table(vfdb_results, label="VFDB hits"),
            ])
        )
    else:
        mo.output.replace(
            mo.md(
                "**VFDB not found.** To enable VFDB enrichment, download the core dataset "
                "protein sequences from https://www.mgc.ac.cn/VFs/download.htm and place "
                f"`VFDB_setA_pro.fas` in `{os.path.dirname(vfdb_fasta)}/`."
            )
        )
    return (vfdb_results,)


@app.cell
def _():
    mo.md("""
    ## Custom Queries

    A pangenome BLAST database is built from the CD-HIT representative
    sequences. Once created, you can search any protein of interest
    against the full pangenome from the command line:

    ```bash
    blastp -query your_query.fa \
           -db data/external/PangenomeDB/PangenomeDB \
           -outfmt 6 -evalue 1e-5
    ```

    Example query sequences can be obtained from
    [UniProt](https://www.uniprot.org/) or
    [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/).
    """)
    return


@app.cell
def _(DATA, SPECIES):
    """Create pangenome BLAST database for custom queries."""
    pangenome_db = os.path.join(DATA, "external/PangenomeDB/PangenomeDB")
    pangenome_input = os.path.join(DATA, "processed/cd-hit-results", SPECIES)

    if os.path.exists(pangenome_input):
        os.makedirs(os.path.dirname(pangenome_db), exist_ok=True)
        if not os.path.exists(pangenome_db + ".phr"):
            make_blast_db(pangenome_input, pangenome_db)
            mo.output.replace(mo.md(f"Created pangenome BLAST DB at `{pangenome_db}`"))
        else:
            mo.output.replace(mo.md(f"Pangenome BLAST DB already exists at `{pangenome_db}`"))
    else:
        mo.output.replace(mo.md(f"Pangenome FASTA not found at `{pangenome_input}`, skipping DB creation"))
    return (pangenome_db,)


@app.cell
def _():
    mo.md("""
    ## Save Results
    """)
    return


@app.cell
def _(OUT, dna_seqs, pangenome_db, prot_seqs, vfdb_results):
    """Save outputs and summary."""
    # Save VFDB results if available
    if vfdb_results is not None:
        vfdb_results.to_csv(os.path.join(OUT, "data/5e_vfdb_blast_results.csv"), index=False)

    # Summary
    summary_rows = [
        {"metric": "representative_protein_seqs", "value": "extracted" if os.path.exists(prot_seqs) else "missing"},
        {"metric": "representative_dna_seqs", "value": "extracted" if os.path.exists(dna_seqs) else "missing"},
        {"metric": "vfdb_hits", "value": str(len(vfdb_results)) if vfdb_results is not None else "skipped"},
        {"metric": "pangenome_db_created", "value": str(os.path.exists(pangenome_db + ".phr"))},
    ]
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUT, "data/5e_blast_summary.csv"), index=False)

    mo.output.replace(
        mo.md(
            "## 5e Summary\n\n"
            + "\n".join(f"- **{r['metric']}**: {r['value']}" for r in summary_rows)
            + f"\n\nSaved to `{os.path.join(OUT, 'data/5e_blast_summary.csv')}`"
        )
    )
    return


if __name__ == "__main__":
    app.run()
