import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml

    from pyphylon.infer_affinities import infer_affinities
    from pyphylon.models import recommended_threshold

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")


@app.cell
def _():
    mo.md("""
    # 5f: Infer Affinities

    Estimate phylon affinities for **new query genomes** that were not part
    of the original NMF decomposition.  The pipeline has two stages:

    1. **MASH cluster classification** — each new strain is compared to all
       reference genomes via MASH distances, and the 5 nearest neighbors
       vote to assign a cluster label.
    2. **NNLS affinity estimation** — Non-Negative Least Squares projects
       each new strain's gene-presence vector onto the NMF basis matrix
       (**L**), producing a phylon-affinity vector that is then binarized
       at 90 % of the original recommended threshold.

    *Prerequisites:* run `workflow/infer_affinities/Snakefile` to generate
    `combined_P_matrix.csv` and `combined_mash_distances.csv`.
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

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    return DATA, FIG, OUT, TEMP


@app.cell
def _():
    mo.md("""
    ## Load Inputs

    Load the NMF basis matrices (**L_norm**, **A_norm**) and the enriched
    metadata table.  The metadata provides the MASH cluster labels used
    for nearest-neighbor classification of new strains.
    """)
    return


@app.cell
def _(DATA, TEMP):
    """Load L_norm, A_norm, and metadata."""
    L_norm = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "L_norm.csv"), index_col=0)
    A_norm = pd.read_csv(os.path.join(DATA, "processed", "nmf-outputs", "A_norm.csv"), index_col=0)

    # Metadata with genome_id as index (needed for MASH cluster lookup)
    metadata = pd.read_csv(
        os.path.join(TEMP, "2d_enriched_metadata.csv"),
        dtype={"genome_id": str},
    )
    metadata = metadata.set_index("genome_id")

    n_phylons = A_norm.shape[0]
    mo.output.replace(
        mo.md(
            f"Loaded NMF inputs:\n\n"
            f"- **L_norm:** {L_norm.shape[0]} genes \u00d7 {L_norm.shape[1]} phylons\n"
            f"- **A_norm:** {A_norm.shape[0]} phylons \u00d7 {A_norm.shape[1]} genomes\n"
            f"- **Metadata:** {len(metadata)} genomes\n"
            f"- **Number of phylons:** {n_phylons}"
        )
    )
    return A_norm, L_norm, metadata, n_phylons


@app.cell
def _():
    mo.md("""
    ## Load Workflow Outputs

    The `workflow/infer_affinities/Snakefile` produces two files:

    - **combined_P_matrix.csv** — gene-presence matrix for new strains,
      built by running CD-HIT-2D against the existing pangenome.
    - **combined_mash_distances.csv** — pairwise MASH distances between
      every reference genome and each new strain.

    If these files are missing, downstream cells skip gracefully.
    """)
    return


@app.cell
def _(DATA):
    """Load infer_affinities workflow outputs (combined P matrix and MASH distances)."""
    infer_dir = os.path.join(DATA, "inferring_affinities")
    p_new_path = os.path.join(infer_dir, "combined_P_matrix.csv")
    mash_dist_path = os.path.join(infer_dir, "combined_mash_distances.csv")

    P_new_raw = None
    mash_distances = None

    if os.path.exists(p_new_path) and os.path.exists(mash_dist_path):
        P_new_raw = pd.read_csv(p_new_path, index_col=0)
        mash_distances = pd.read_csv(mash_dist_path, dtype="object").set_index("genome_id").astype(float)
        mo.output.replace(
            mo.md(
                f"Loaded workflow outputs:\n\n"
                f"- **P_new:** {P_new_raw.shape[0]} genes \u00d7 {P_new_raw.shape[1]} new genomes\n"
                f"- **MASH distances:** {mash_distances.shape[0]} ref genomes \u00d7 {mash_distances.shape[1]} new genomes\n"
                f"- **New strains:** {list(mash_distances.columns)}"
            )
        )
    else:
        missing = []
        if not os.path.exists(p_new_path):
            missing.append(f"`{p_new_path}`")
        if not os.path.exists(mash_dist_path):
            missing.append(f"`{mash_dist_path}`")
        mo.output.replace(
            mo.md(
                f"**Workflow outputs not found** \u2014 skipping affinity inference.\n\n"
                f"Missing: {', '.join(missing)}\n\n"
                f"Run `workflow/infer_affinities/Snakefile` first to generate these files."
            )
        )
    return P_new_raw, mash_distances


@app.cell
def _():
    mo.md("""
    ## MASH Cluster Assignments

    Each new strain is classified into a MASH cluster using **5-nearest-neighbor
    majority vote**: the five reference genomes with the smallest MASH distance
    are identified, and the most common cluster label among them is assigned.
    """)
    return


@app.cell
def _(mash_distances, metadata):
    """Classify new strains into MASH clusters via 5-NN majority vote."""
    strain_mash_clusters = None

    if mash_distances is not None:
        strain_mash_clusters = pd.DataFrame(columns=["best_cluster", "nearest_mash_distance_mean"])

        for new_strain in mash_distances.columns:
            top_distances = mash_distances[new_strain].sort_values().index[:5]
            mash_cluster = int(float(metadata.loc[top_distances].mash_cluster.value_counts().idxmax()))
            strain_mash_clusters.loc[new_strain] = [
                mash_cluster,
                mash_distances[new_strain].sort_values().head(5).mean(),
            ]

        mo.output.replace(
            mo.ui.table(strain_mash_clusters.reset_index().rename(columns={"index": "strain"}))
        )
    else:
        mo.output.replace(mo.md("**Skipped** \u2014 workflow outputs not available."))
    return (strain_mash_clusters,)


@app.cell
def _():
    mo.md("""
    ## Infer Phylon Affinities

    Phylon affinities are estimated via **Non-Negative Least Squares (NNLS)**:
    each new strain's gene-presence vector is projected onto the NMF basis
    matrix **L_norm**, yielding a continuous affinity score per phylon.

    The resulting matrix is binarized at **90 %** of the `recommended_threshold`
    from the original **A_norm** matrix.  The 0.9 multiplier compensates for
    CD-HIT clustering incongruities when mapping new strains against an
    existing pangenome.
    """)
    return


@app.cell
def _(A_norm, L_norm, P_new_raw, mash_distances):
    """Infer affinities via NNLS and binarize."""
    A_new = None
    A_new_binarized = None

    if P_new_raw is not None and mash_distances is not None:
        # Gene alignment: reindex P_new to L_norm's gene set
        P_new_aligned = P_new_raw.reindex(L_norm.index, fill_value=0)

        # Infer affinities via NNLS
        A_new_arr = infer_affinities(L_norm.to_numpy(), P_new_aligned.to_numpy(), n_jobs=4)
        A_new = pd.DataFrame(A_new_arr, index=L_norm.columns, columns=P_new_aligned.columns)

        # Binarize using 90% of recommended threshold from original A_norm
        A_new_binarized = pd.DataFrame(np.zeros_like(A_new.values), index=A_new.index, columns=A_new.columns)
        for idx in A_new.index:
            phylon_num = idx.replace("phylon", "")
            threshold = recommended_threshold(A_norm, phylon_num) * 0.9
            cond = A_new.loc[idx] >= threshold
            A_new_binarized.loc[idx, cond] = 1

        mo.output.replace(
            mo.md(
                f"Affinity inference results:\n\n"
                f"- **New strains:** {A_new.shape[1]}\n"
                f"- **Phylons:** {A_new.shape[0]}\n"
                f"- **Active affinities (binarized):** {int(A_new_binarized.sum().sum())} "
                f"across {int((A_new_binarized.sum(axis=0) > 0).sum())} strains"
            )
        )
    else:
        mo.output.replace(mo.md("**Skipped** \u2014 workflow outputs not available."))
    return A_new, A_new_binarized


@app.cell
def _():
    mo.md("""
    ## Affinity Heatmap

    Continuous affinity scores (greyscale) for each new strain across all
    phylons.  Darker cells indicate stronger affinity.
    """)
    return


@app.cell
def _(A_new, FIG):
    """Generate affinity heatmap."""
    if A_new is not None:
        fig_hm, ax_hm = plt.subplots(
            figsize=(max(6, A_new.shape[1] * 0.8), max(4, A_new.shape[0] * 0.3))
        )
        sns.heatmap(A_new, cmap="Greys", ax=ax_hm)
        ax_hm.set_title("Inferred Affinities (New Strains)")
        ax_hm.set_ylabel("Phylon")
        ax_hm.set_xlabel("Genome")
        plt.tight_layout()
        fig_hm.savefig(os.path.join(FIG, "5f_affinity_heatmap.png"), bbox_inches="tight")
        mo.output.replace(fig_hm)
    else:
        mo.output.replace(mo.md("**Skipped** \u2014 no affinity data to visualize."))
    return


@app.cell
def _():
    mo.md("""
    ## Save Results
    """)
    return


@app.cell
def _(A_new, A_new_binarized, DATA, OUT, n_phylons, strain_mash_clusters):
    """Save affinity matrices, MASH clusters, and summary CSV."""
    nmf_dir = os.path.join(DATA, "processed", "nmf-outputs")
    summary_rows = []

    if A_new is not None:
        A_new.to_csv(os.path.join(nmf_dir, "A_new.csv"))
        A_new_binarized.to_csv(os.path.join(nmf_dir, "A_new_binarized.csv"))
        strain_mash_clusters.to_csv(os.path.join(OUT, "data", "5f_mash_clusters.csv"))

        summary_rows = [
            {"metric": "n_new_strains", "value": str(A_new.shape[1])},
            {"metric": "n_phylons", "value": str(n_phylons)},
            {"metric": "active_affinities", "value": str(int(A_new_binarized.sum().sum()))},
            {"metric": "strains_with_affinity", "value": str(int((A_new_binarized.sum(axis=0) > 0).sum()))},
            {"metric": "mean_affinities_per_strain", "value": f"{A_new_binarized.sum(axis=0).mean():.2f}"},
            {"metric": "status", "value": "completed"},
        ]
    else:
        summary_rows = [
            {"metric": "status", "value": "skipped"},
            {"metric": "reason", "value": "workflow outputs not found"},
        ]

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUT, "data", "5f_affinity_summary.csv"), index=False)

    mo.output.replace(
        mo.md("### Summary\n\n" + "\n".join(f"- **{r['metric']}**: {r['value']}" for r in summary_rows))
    )
    return


if __name__ == "__main__":
    app.run()
