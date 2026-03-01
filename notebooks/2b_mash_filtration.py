import marimo

__generated_with = "0.20.2"
app = marimo.App(width="medium")

with app.setup:
    import os
    import sys

    import marimo as mo
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    import yaml

    from pyphylon.mash import cluster_corr_dist, remove_bad_strains, sensitivity_analysis


@app.cell
def _():
    """Parse config and set up directories."""
    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    TEMP = CONFIG["TEMP_DIR"]
    OUTPUT = CONFIG["OUTPUT_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    SMALL_CLUSTER_LIMIT = CONFIG.get("SMALL_CLUSTER_LIMIT", 5)

    os.makedirs(TEMP, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT, "data"), exist_ok=True)

    return CONFIG, FIG, OUTPUT, SMALL_CLUSTER_LIMIT, TEMP


@app.cell
def _(TEMP):
    """Load 2a outputs: genome summary and metadata."""
    input_summary = pd.read_csv(
        os.path.join(TEMP, "2a_genome_summary.csv"),
        dtype={"genome_id": str},
    )
    input_metadata = pd.read_csv(
        os.path.join(TEMP, "2a_genome_metadata.csv"),
        dtype={"genome_id": str},
    )

    mo.output.replace(
        mo.md(
            f"Loaded from 2a:\n\n"
            f"- **genome_summary:** {input_summary.shape[0]} rows\n"
            f"- **genome_metadata:** {input_metadata.shape[0]} rows"
        )
    )
    return input_metadata, input_summary


@app.cell
def _(TEMP):
    """Load mash_distances.txt and pivot to square distance matrix."""
    mash_dist_path = os.path.join(TEMP, "2b_mash", "mash_distances.txt")
    names = ["genome1", "genome2", "mash_distance", "p_value", "matching_hashes"]

    df_mash = pd.read_csv(mash_dist_path, sep="\t", names=names)
    # Extract genome IDs from file paths (e.g. "temp/.../1201032.3.fna" -> "1201032.3")
    df_mash["genome1"] = df_mash["genome1"].apply(lambda x: x.split("/")[-1].split(".fna")[0])
    df_mash["genome2"] = df_mash["genome2"].apply(lambda x: x.split("/")[-1].split(".fna")[0])

    mash_square_raw = df_mash.pivot(index="genome1", columns="genome2", values="mash_distance")

    mo.output.replace(mo.md(f"Mash square matrix: **{mash_square_raw.shape[0]}** x **{mash_square_raw.shape[1]}**"))
    return (mash_square_raw,)


@app.cell
def _(FIG, mash_square_raw):
    """Plot raw mash distance heatmap."""
    fig_heatmap, ax_heatmap = plt.subplots(figsize=(8, 6))
    sns.heatmap(mash_square_raw, cmap="viridis", ax=ax_heatmap, xticklabels=False, yticklabels=False)
    ax_heatmap.set_title("Raw Mash Distance Heatmap")
    fig_heatmap.savefig(os.path.join(FIG, "2b_mash_heatmap.png"), bbox_inches="tight", dpi=200)
    mo.output.replace(fig_heatmap)
    return


@app.cell
def _(FIG, input_metadata, mash_square_raw):
    """Compute correlation distance matrix (cached) and apply filtration steps."""
    # --- Step 1: Pearson correlation distance (EXPENSIVE, cached) ---
    with mo.persistent_cache("2b_corr_matrix"):
        mash_corr_dist_raw = 1 - mash_square_raw.corr()

    # --- Step 2: Filter to 2a (scrubbed) strains ---
    scrubbed_strains = input_metadata["genome_id"].astype(str)
    scrubbed_strains = scrubbed_strains[scrubbed_strains.isin(mash_square_raw.index.astype(str))]

    mash_square_scrubbed = mash_square_raw.loc[scrubbed_strains, scrubbed_strains]
    mash_corr_dist_scrubbed = mash_corr_dist_raw.loc[scrubbed_strains, scrubbed_strains]

    # --- Step 3: Auto-detect reference/representative strains ---
    ref_mask = input_metadata["reference_genome"].notna()
    ref_ids = input_metadata.loc[ref_mask, "genome_id"].astype(str)
    repr_strains = sorted(ref_ids[ref_ids.isin(mash_square_scrubbed.index)].tolist())

    if not repr_strains:
        # Fallback: medoid (genome with smallest mean mash distance)
        repr_strains = [mash_square_scrubbed.mean(axis=1).idxmin()]

    # --- Step 4: Filter by mash distance < 0.05 to reference strains ---
    current_square = mash_square_scrubbed
    current_corr_dist = mash_corr_dist_scrubbed
    for repr_strain in repr_strains:
        cond = current_square.loc[repr_strain] < 0.05
        good_strains = current_square.loc[repr_strain][cond].index
        current_square = current_square.loc[good_strains, good_strains]
        current_corr_dist = current_corr_dist.loc[good_strains, good_strains]

    mash_square_filtered = current_square
    mash_corr_dist_filtered = current_corr_dist

    # --- Histogram: pre-filtration ---
    fig_pre, ax_pre = plt.subplots()
    sns.histplot(mash_square_scrubbed.values.flatten(), ax=ax_pre)
    ax_pre.set_title(f"All mash distances ({mash_square_scrubbed.shape[0]} genomes)")
    fig_pre.savefig(os.path.join(FIG, "2b_mash_dist_all.png"), bbox_inches="tight", dpi=200)

    # --- Histogram: post-filtration ---
    fig_post, ax_post = plt.subplots()
    sns.histplot(mash_square_filtered.values.flatten(), ax=ax_post)
    ax_post.set_title(f"Post-filtration mash distances ({mash_square_filtered.shape[0]} genomes)")
    fig_post.savefig(os.path.join(FIG, "2b_mash_dist_filtered.png"), bbox_inches="tight", dpi=200)

    mo.output.replace(
        mo.vstack(
            [
                mo.md(
                    f"Correlation matrix: **{mash_corr_dist_raw.shape[0]}** x **{mash_corr_dist_raw.shape[1]}**\n\n"
                    f"After 2a strain filter: **{mash_square_scrubbed.shape[0]}** genomes\n\n"
                    f"Representative strains: {repr_strains}\n\n"
                    f"After mash < 0.05 filter: **{mash_square_filtered.shape[0]}** genomes"
                ),
                fig_pre,
                fig_post,
            ]
        )
    )
    return mash_corr_dist_filtered, mash_square_filtered


@app.cell
def _(FIG, input_summary, mash_corr_dist_filtered, mash_square_filtered):
    """Extract complete genomes and run sensitivity analysis for clustering threshold."""
    # Subset to complete genomes
    cond_complete = input_summary["genome_status"] == "Complete"
    complete_ids = set(input_summary[cond_complete]["genome_id"].astype(str))
    complete_seqs = sorted(complete_ids.intersection(set(mash_square_filtered.index)))

    mash_square_complete = mash_square_filtered.loc[complete_seqs, complete_seqs]
    mash_corr_dist_complete = mash_corr_dist_filtered.loc[complete_seqs, complete_seqs]

    # Sensitivity analysis: find elbow threshold
    tmp_sa, df_temp_sa, elbow_idx, elbow_threshold = sensitivity_analysis(mash_corr_dist_complete)

    # Plot
    fig_sa, ax_sa = plt.subplots(figsize=(5, 3.5))
    ax_sa.plot(tmp_sa["threshold"], tmp_sa["num_clusters"])
    ax_sa.axhline(y=df_temp_sa["num_clusters"][elbow_idx], c="#ff00ff", linestyle="--")
    ax_sa.set_ylabel("num_clusters")
    ax_sa.set_xlabel("threshold")
    fig_sa.suptitle(
        f"Elbow at {df_temp_sa['num_clusters'][elbow_idx]} clusters (threshold: {elbow_threshold:.4f})",
        y=1.02,
    )
    fig_sa.tight_layout()
    fig_sa.savefig(os.path.join(FIG, "2b_sensitivity_analysis.png"), bbox_inches="tight", dpi=200)

    # Round up threshold (following original notebook convention)
    elbow_threshold = elbow_threshold + 0.1

    mo.output.replace(
        mo.vstack(
            [
                mo.md(
                    f"Complete genomes in filtered set: **{len(complete_seqs)}**\n\n"
                    f"Elbow threshold (+ 0.1 offset): **{elbow_threshold:.4f}**"
                ),
                fig_sa,
            ]
        )
    )
    return elbow_threshold, mash_corr_dist_complete, mash_square_complete


@app.cell
def _(FIG, elbow_threshold, mash_corr_dist_complete, mash_square_complete):
    """Initial clustering and clustermap visualization."""
    initial_link, _initial_dist, initial_clst = cluster_corr_dist(mash_corr_dist_complete, thresh=elbow_threshold)

    # Color each cluster
    cm = matplotlib.colormaps.get_cmap("tab20")
    initial_clr = dict(zip(sorted(initial_clst["cluster"].unique()), cm.colors, strict=False))
    initial_clst["color"] = initial_clst["cluster"].map(initial_clr)

    # Clustermap
    legend_patches = [mpatches.Patch(color=c, label=str(k)) for k, c in initial_clr.items()]
    sns.set(rc={"figure.facecolor": "white"})
    g_initial = sns.clustermap(
        mash_square_complete,
        figsize=(6, 6),
        row_linkage=initial_link,
        col_linkage=initial_link,
        col_colors=initial_clst["color"],
        yticklabels=False,
        xticklabels=False,
        cmap="BrBG_r",
        robust=True,
    )
    g_initial.ax_heatmap.legend(
        loc="upper left",
        bbox_to_anchor=(1.01, 0.85),
        handles=legend_patches,
        frameon=True,
        title="Clusters",
    )
    g_initial.savefig(os.path.join(FIG, "2b_clustermap_initial.png"), bbox_inches="tight", dpi=200)

    # Cluster size bar chart
    fig_sizes, ax_sizes = plt.subplots()
    initial_clst["cluster"].value_counts().sort_index().plot.bar(ax=ax_sizes)
    ax_sizes.set_xlabel("Cluster")
    ax_sizes.set_ylabel("Count")
    ax_sizes.set_title("Initial cluster sizes")
    fig_sizes.tight_layout()
    fig_sizes.savefig(os.path.join(FIG, "2b_cluster_sizes_initial.png"), bbox_inches="tight", dpi=200)

    mo.output.replace(
        mo.vstack(
            [
                mo.md(
                    f"Initial clustering: **{len(initial_clr)}** clusters\n\n"
                    f"Smallest cluster: **{initial_clst['cluster'].value_counts().min()}** genomes"
                ),
                g_initial.figure,
                fig_sizes,
            ]
        )
    )
    return


@app.cell
def _(SMALL_CLUSTER_LIMIT, elbow_threshold, mash_corr_dist_complete, mash_square_complete):
    """Iteratively remove small clusters until all clusters >= SMALL_CLUSTER_LIMIT."""
    # Work on copies to avoid mutating upstream cell outputs
    iter_square = mash_square_complete.copy()
    iter_corr_dist = mash_corr_dist_complete.copy()

    # Initial clustering to seed the loop
    _link, _dist, iter_clst = cluster_corr_dist(iter_corr_dist, thresh=elbow_threshold)

    iteration = 0
    prev_n_clusters = 0
    curr_n_clusters = len(iter_clst["cluster"].unique())
    log_lines = []

    while abs(prev_n_clusters - curr_n_clusters) > 0:
        iteration += 1
        log_lines.append(f"Iteration {iteration}: {curr_n_clusters} clusters, {iter_square.shape[0]} genomes")

        # Identify small clusters
        cluster_counts = iter_clst["cluster"].value_counts()
        bad_cluster_ids = cluster_counts[cluster_counts < SMALL_CLUSTER_LIMIT].index

        # Identify genomes in small clusters
        bad_genomes = [genome for genome in iter_square.index if iter_clst.loc[genome, "cluster"] in bad_cluster_ids]

        # Remove bad genomes from both matrices
        iter_square = remove_bad_strains(iter_square, bad_genomes)
        iter_corr_dist = remove_bad_strains(iter_corr_dist, bad_genomes)

        # Re-cluster
        prev_n_clusters = curr_n_clusters
        _link, _dist, iter_clst = cluster_corr_dist(iter_corr_dist, thresh=elbow_threshold)
        curr_n_clusters = len(iter_clst["cluster"].unique())

    log_lines.append(
        f"Converged after {iteration} iterations: {curr_n_clusters} clusters, {iter_square.shape[0]} genomes"
    )

    # Final clustering to get definitive linkage + cluster assignments
    final_link, _final_dist, final_clst = cluster_corr_dist(iter_corr_dist, thresh=elbow_threshold)

    # Assert minimum cluster size
    min_cluster_size = final_clst["cluster"].value_counts().min()
    assert min_cluster_size >= SMALL_CLUSTER_LIMIT, f"Minimum cluster size {min_cluster_size} < {SMALL_CLUSTER_LIMIT}"

    mash_square_final = iter_square
    mash_corr_dist_final = iter_corr_dist

    mo.output.replace(
        mo.md(
            "### Iterative small-cluster removal\n\n"
            + "\n".join(f"- {line}" for line in log_lines)
            + f"\n\nMin cluster size: **{min_cluster_size}**"
        )
    )
    return final_clst, final_link, mash_corr_dist_final, mash_square_final


@app.cell
def _(FIG, final_clst, final_link, mash_square_final):
    """Final clustermap and cluster size visualization."""
    # Color clusters
    cm_final = matplotlib.colormaps.get_cmap("tab20")
    final_clr = dict(zip(sorted(final_clst["cluster"].unique()), cm_final.colors, strict=False))
    final_clst_colored = final_clst.copy()
    final_clst_colored["color"] = final_clst_colored["cluster"].map(final_clr)

    # Clustermap
    legend_patches_final = [mpatches.Patch(color=c, label=str(k)) for k, c in final_clr.items()]
    sns.set(rc={"figure.facecolor": "white"})
    g_final = sns.clustermap(
        mash_square_final,
        figsize=(6, 6),
        row_linkage=final_link,
        col_linkage=final_link,
        col_colors=final_clst_colored["color"],
        yticklabels=False,
        xticklabels=False,
        cmap="BrBG_r",
        robust=True,
    )
    g_final.ax_heatmap.legend(
        loc="upper left",
        bbox_to_anchor=(1.05, 0.85),
        handles=legend_patches_final,
        frameon=True,
        title="Clusters",
    )
    g_final.savefig(os.path.join(FIG, "2b_clustermap_final.png"), bbox_inches="tight", dpi=200)

    # Cluster size bar chart
    fig_final_sizes, ax_final_sizes = plt.subplots()
    final_clst["cluster"].value_counts().sort_index().plot.bar(ax=ax_final_sizes)
    ax_final_sizes.set_xlabel("Cluster")
    ax_final_sizes.set_ylabel("Count")
    ax_final_sizes.set_title("Final cluster sizes")
    fig_final_sizes.tight_layout()
    fig_final_sizes.savefig(os.path.join(FIG, "2b_cluster_sizes_final.png"), bbox_inches="tight", dpi=200)

    mo.output.replace(mo.vstack([g_final.figure, fig_final_sizes]))
    return


@app.cell
def _(final_clst, input_metadata, input_summary):
    """Filter metadata/summary to surviving genomes and add mash_cluster column."""
    survived = final_clst.index
    final_metadata = input_metadata[input_metadata["genome_id"].astype(str).isin(survived)].reset_index(drop=True)
    final_summary = input_summary[input_summary["genome_id"].astype(str).isin(survived)].reset_index(drop=True)

    # Add mash_cluster column
    cluster_map = final_clst["cluster"].to_dict()
    final_metadata["mash_cluster"] = final_metadata["genome_id"].astype(str).map(cluster_map)

    mo.output.replace(
        mo.md(
            f"Final dataset:\n\n"
            f"- **Genomes:** {final_metadata.shape[0]}\n"
            f"- **Clusters:** {final_metadata['mash_cluster'].nunique()}\n"
            f"- **Summary rows:** {final_summary.shape[0]}"
        )
    )
    return final_metadata, final_summary


@app.cell
def _(OUTPUT, TEMP, final_metadata, final_summary, mash_corr_dist_filtered, mash_square_filtered):
    """Save all CSV outputs."""
    # Temp outputs (consumed by downstream notebooks)
    final_summary.to_csv(os.path.join(TEMP, "2b_genome_summary.csv"), index=False)
    final_metadata.to_csv(os.path.join(TEMP, "2b_genome_metadata.csv"), index=False)

    # Permanent outputs (square matrices — keep index)
    mash_square_filtered.to_csv(os.path.join(OUTPUT, "data", "2b_mash_square.csv"))
    mash_corr_dist_filtered.to_csv(os.path.join(OUTPUT, "data", "2b_mash_corr_dist.csv"))

    mo.output.replace(
        mo.md(
            f"Saved:\n\n"
            f"- `{TEMP}2b_genome_summary.csv` ({final_summary.shape[0]} rows)\n"
            f"- `{TEMP}2b_genome_metadata.csv` ({final_metadata.shape[0]} rows)\n"
            f"- `{OUTPUT}data/2b_mash_square.csv` ({mash_square_filtered.shape[0]} x {mash_square_filtered.shape[1]})\n"
            f"- `{OUTPUT}data/2b_mash_corr_dist.csv` ({mash_corr_dist_filtered.shape[0]} x {mash_corr_dist_filtered.shape[1]})"
        )
    )
    return


if __name__ == "__main__":
    app.run()
