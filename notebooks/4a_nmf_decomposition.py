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
    import prince
    import scipy.spatial.distance
    import seaborn as sns
    import yaml
    from scipy.cluster import hierarchy as hc
    from scipy.cluster.hierarchy import cophenet

    from pyphylon.models import (
        binarize_nmf_outputs,
        calculate_nmf_reconstruction_metrics,
        generate_nmf_reconstructions,
        normalize_nmf_outputs,
        run_nmf,
    )
    from pyphylon.pangenome import connectivity, get_gene_frequency_submatrices

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["figure.dpi"] = 200
    sns.set_palette("deep")
    sns.set_context("paper")
    sns.set_style("whitegrid")

    def build_consensus(P_submatrices, best_A_binarized_dict, genome_ids, exclude_pairs=None, thresh=0.7):
        """Build consensus matrix, cluster, and compute cophenetic correlation.

        Returns (df_consensus, clusters, linkage, coph_cor).
        """
        exclude_pairs = exclude_pairs or set()
        n = len(genome_ids)
        consensus = np.zeros((n, n))
        count = 0
        for min_key in sorted(P_submatrices.keys()):
            for max_key in sorted(P_submatrices.keys()):
                if min_key >= max_key:
                    continue
                if P_submatrices[min_key].get(max_key) is None:
                    continue
                if (min_key, max_key) in exclude_pairs:
                    continue
                conn = connectivity(
                    P_submatrices[min_key][max_key].values,
                    best_A_binarized_dict[(min_key, max_key)].values,
                )
                consensus += conn
                count += 1
        consensus /= count

        df_consensus = pd.DataFrame(consensus, index=genome_ids, columns=genome_ids)
        df_dist = 1 - df_consensus
        link = hc.linkage(scipy.spatial.distance.squareform(df_dist), method="ward")
        dist_sq = scipy.spatial.distance.squareform(df_dist)
        clusters = pd.DataFrame(index=genome_ids)
        clusters["cluster"] = hc.fcluster(link, thresh * dist_sq.max(), "distance")

        # Cophenetic correlation
        avec = np.array([consensus[i, j] for i in range(n - 1) for j in range(i + 1, n)])
        Y = 1 - avec
        Z = hc.linkage(Y, method="ward")
        coph_cor, _ = cophenet(Z, Y)

        return df_consensus, clusters, link, coph_cor


@app.cell
def _():
    mo.md(
        """
        # 4a: NMF Decomposition — Phylon Identification

        Non-negative Matrix Factorization (NMF) decomposes the accessory genome
        presence/absence matrix **P** into two lower-rank matrices **L** (genomes × phylons)
        and **A** (phylons × genes), identifying co-occurring gene sets (phylons).

        **Workflow:**

        1. **MCA rank selection** — Multiple Correspondence Analysis identifies candidate
           ranks by explained-variance thresholds (70–90 %).
        2. **NMF screening** — Initial decomposition at MCA-derived ranks, then refinement
           around the best AIC rank.
        3. **Submatrix NMF** — Repeat decomposition on gene-frequency submatrices (e.g.
           0–25 %, 25–50 %, …) to assess robustness.
        4. **Consensus clustering** — Average connectivity matrices across submatrices and
           cluster with Ward's method; cophenetic correlation measures cluster stability.
        """
    )


@app.cell
def _():
    mo.md(
        """
        ## Setup
        """
    )


@app.cell
def _():
    """Parse config and set up directories."""
    config_path = "config.yml"
    if "--config" in sys.argv:
        config_path = sys.argv[sys.argv.index("--config") + 1]

    with open(config_path) as f:
        CONFIG = yaml.safe_load(f)

    TEMP = CONFIG["TEMP_DIR"]
    DATA = CONFIG["DATA_DIR"]
    FIG = CONFIG["FIGURES_DIR"]
    OUT = CONFIG["OUTPUT_DIR"]
    SPECIES = CONFIG["PG_NAME"]
    MASH_RANK = CONFIG.get("MASH_RANK", 16)

    os.makedirs(FIG, exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    os.makedirs(os.path.join(DATA, "processed", "nmf-outputs"), exist_ok=True)

    return CONFIG, DATA, FIG, MASH_RANK, OUT, SPECIES, TEMP


@app.cell
def _(DATA, TEMP):
    """Load accessory genome matrix, submatrices, and metadata."""
    df_acc = pd.read_csv(os.path.join(DATA, "processed", "CAR_genomes", "df_acc.csv"), index_col=0)

    P_submatrices = get_gene_frequency_submatrices(df_acc)
    P_submatrices[0][100] = None  # remove full accessory (redundant with df_acc)

    metadata = pd.read_csv(
        os.path.join(TEMP, "2d_enriched_metadata.csv"),
        dtype={"genome_id": str},
    )

    mo.output.replace(
        mo.md(
            f"Loaded inputs:\n\n"
            f"- **Accessory genome:** {df_acc.shape[0]} genes x {df_acc.shape[1]} genomes\n"
            f"- **Submatrices:** {sum(1 for mk in P_submatrices for xk in P_submatrices[mk] if mk < xk and P_submatrices[mk].get(xk) is not None)} frequency-band pairs\n"
            f"- **Metadata:** {metadata.shape[0]} genomes"
        )
    )
    return df_acc, P_submatrices, metadata


@app.cell
def _(MASH_RANK, df_acc):
    """MCA for rank analysis (cached — ~5 min)."""
    with mo.persistent_cache("4a_mca"):
        mca = prince.MCA(
            n_components=df_acc.shape[1],
            n_iter=3,
            copy=True,
            check_input=True,
            engine="sklearn",
            random_state=42,
        ).fit(df_acc)
        cumulative_variance = pd.Series(mca.percentage_of_variance_).cumsum()
        threshold = {n: int(cumulative_variance[cumulative_variance >= n].index[0]) for n in range(1, 99)}
        rank_list = sorted(
            set(
                [
                    2,
                    MASH_RANK,
                    threshold[70],
                    threshold[75],
                    threshold[80],
                    threshold[85],
                    threshold[90],
                ]
            )
        )

    mo.output.replace(mo.md(f"MCA complete. Rank candidates: **{rank_list}**"))
    return cumulative_variance, rank_list, threshold


@app.cell
def _():
    mo.md(
        """
        ## MCA Cumulative Variance

        MCA captures the dominant axes of variation in the binary presence/absence
        matrix. Vertical lines mark the number of components needed to explain
        70 %, 75 %, 80 %, 85 %, and 90 % of the total variance — these become
        candidate ranks for NMF.
        """
    )


@app.cell
def _(FIG, cumulative_variance, threshold):
    """Plot MCA cumulative variance with rank thresholds."""
    n_significant = (cumulative_variance.diff().fillna(cumulative_variance.iloc[0]) > 0.01).sum()

    fig_mca, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(16, 5))

    # Full variance plot
    ax_full.plot(cumulative_variance, marker=".", markersize=2, color="blue")
    ax_full.set_title("Cumulative Explained Variance (Full)")
    ax_full.set_xlabel("Number of Dimensions")
    ax_full.set_ylabel("Cumulative Explained Variance (%)")

    # Zoomed plot (significant components only)
    ax_zoom.plot(cumulative_variance[:n_significant], marker=".", markersize=2, color="blue")
    ax_zoom.set_title(f"Cumulative Explained Variance (First {n_significant} Components)")
    ax_zoom.set_xlabel("Number of Dimensions")
    ax_zoom.set_ylabel("Cumulative Explained Variance (%)")

    # Add threshold lines to both axes
    colors = {"70": "grey", "75": "limegreen", "80": "purple", "85": "pink", "90": "maroon"}
    for pct_str, color in colors.items():
        pct = int(pct_str)
        for ax in [ax_full, ax_zoom]:
            ax.axvline(x=threshold[pct], color=color, linestyle="--", label=f"{pct}%: {threshold[pct]}")

    ax_full.legend(fontsize=8)
    ax_zoom.legend(fontsize=8)
    plt.tight_layout()

    fig_mca.savefig(os.path.join(FIG, "4a_mca_variance.png"), bbox_inches="tight")
    mo.output.replace(fig_mca)


@app.cell
def _():
    mo.md(
        """
        ## NMF Decomposition

        NMF is run in two phases:

        1. **Initial screening** across MCA-derived ranks to find the best AIC.
        2. **Refinement** — re-run with additional ranks near the best AIC model.

        The same procedure is then repeated on each gene-frequency submatrix
        (e.g. 0–25 %, 25–50 %) to collect per-submatrix best models for
        consensus clustering.
        """
    )


@app.cell
def _(df_acc, rank_list):
    """Full-genome NMF: initial screening + refinement around best AIC rank (cached — ~30 min)."""
    with mo.persistent_cache("4a_nmf_full"):
        # Phase 1: Initial screening across MCA-derived ranks
        W_dict, H_dict = run_nmf(data=df_acc, ranks=rank_list, max_iter=50_000)
        L_norm_dict, A_norm_dict = normalize_nmf_outputs(df_acc, W_dict, H_dict)
        L_binarized_dict, A_binarized_dict = binarize_nmf_outputs(L_norm_dict, A_norm_dict)
        P_recon, P_err, P_conf = generate_nmf_reconstructions(df_acc, L_binarized_dict, A_binarized_dict)
        df_metrics = calculate_nmf_reconstruction_metrics(P_recon, P_conf)

        # Phase 2: Refine around best AIC rank
        best_rank_initial = df_metrics["AIC"].idxmin()
        extra_ranks = sorted(set(rank_list + list(range(max(2, best_rank_initial - 3), best_rank_initial + 4))))
        W_dict, H_dict = run_nmf(data=df_acc, ranks=extra_ranks, max_iter=50_000)
        L_norm_dict, A_norm_dict = normalize_nmf_outputs(df_acc, W_dict, H_dict)
        L_binarized_dict, A_binarized_dict = binarize_nmf_outputs(L_norm_dict, A_norm_dict)
        P_recon, P_err, P_conf = generate_nmf_reconstructions(df_acc, L_binarized_dict, A_binarized_dict)
        df_metrics = calculate_nmf_reconstruction_metrics(P_recon, P_conf)

    _best_rank = df_metrics["AIC"].idxmin()
    _metrics_sorted = df_metrics.sort_values("AIC")
    _metrics_header = "| Rank | " + " | ".join(_metrics_sorted.columns) + " |\n"
    _metrics_sep = "|---" * (len(_metrics_sorted.columns) + 1) + "|\n"
    _metrics_rows = "\n".join(
        f"| {idx} | " + " | ".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in row) + " |"
        for idx, row in _metrics_sorted.iterrows()
    )
    mo.output.replace(
        mo.md(
            f"NMF complete. Best rank by AIC: **{_best_rank}** (from {len(extra_ranks)} candidates)\n\n"
            + _metrics_header
            + _metrics_sep
            + _metrics_rows
        )
    )
    return A_norm_dict, A_binarized_dict, L_norm_dict, L_binarized_dict, df_metrics, extra_ranks


@app.cell
def _(A_binarized_dict, A_norm_dict, L_binarized_dict, L_norm_dict, P_submatrices, df_acc, df_metrics, extra_ranks):
    """Submatrix NMF decomposition + best model extraction (cached — ~10-15 hrs)."""
    with mo.persistent_cache("4a_nmf_submatrix"):
        # Initialize with best full-genome model
        best_rank_full = df_metrics["AIC"].idxmin()
        best_ranks_dict = {"full": best_rank_full}
        best_L_norm_dict = {"full": L_norm_dict[best_rank_full]}
        best_A_norm_dict = {"full": A_norm_dict[best_rank_full]}
        best_L_binarized_dict = {"full": L_binarized_dict[best_rank_full]}
        best_A_binarized_dict = {"full": A_binarized_dict[best_rank_full]}

        for min_key in sorted(P_submatrices.keys()):
            for max_key in sorted(P_submatrices.keys()):
                if min_key >= max_key:
                    continue
                if P_submatrices[min_key].get(max_key) is None:
                    continue
                sub = P_submatrices[min_key][max_key]
                if sub.empty:
                    continue

                # Filter ranks to valid range for this submatrix
                valid_ranks = [r for r in extra_ranks if r < min(sub.shape)]
                if not valid_ranks:
                    continue

                W_sub, H_sub = run_nmf(data=sub, ranks=valid_ranks, max_iter=50_000)
                L_norm_sub, A_norm_sub = normalize_nmf_outputs(sub, W_sub, H_sub)
                L_bin_sub, A_bin_sub = binarize_nmf_outputs(L_norm_sub, A_norm_sub)
                P_recon_sub, _, P_conf_sub = generate_nmf_reconstructions(sub, L_bin_sub, A_bin_sub)
                df_met_sub = calculate_nmf_reconstruction_metrics(P_recon_sub, P_conf_sub)

                sub_best = df_met_sub["AIC"].idxmin()
                best_ranks_dict[(min_key, max_key)] = sub_best
                best_L_norm_dict[(min_key, max_key)] = L_norm_sub[sub_best]
                best_A_norm_dict[(min_key, max_key)] = A_norm_sub[sub_best]
                best_L_binarized_dict[(min_key, max_key)] = L_bin_sub[sub_best]
                best_A_binarized_dict[(min_key, max_key)] = A_bin_sub[sub_best]

    mo.output.replace(
        mo.md(
            "Submatrix NMF complete. Best ranks per submatrix:\n\n"
            + "\n".join(f"- **{k}**: rank {v}" for k, v in sorted(best_ranks_dict.items(), key=str))
        )
    )
    return best_ranks_dict, best_L_norm_dict, best_A_norm_dict, best_L_binarized_dict, best_A_binarized_dict


@app.cell
def _():
    mo.md(
        """
        ## Consensus Matrix

        Connectivity matrices from each submatrix decomposition are averaged to
        build a consensus matrix. Ward's-method clustering is applied, and the
        cophenetic correlation coefficient measures how faithfully the dendrogram
        preserves pairwise distances (ideally ≥ 0.7).
        """
    )


@app.cell
def _(P_submatrices, best_A_binarized_dict, df_acc):
    """Build consensus matrices (full and filtered) with cophenetic correlation."""
    # Full consensus (exclude only (0,100) which is None)
    df_consensus, consensus_clst, link, coph_cor = build_consensus(P_submatrices, best_A_binarized_dict, df_acc.columns)

    # Filtered consensus (additionally exclude (50,100) and (75,100))
    df_consensus_filt, consensus_clst_filt, link_filt, coph_cor_filt = build_consensus(
        P_submatrices, best_A_binarized_dict, df_acc.columns, exclude_pairs={(50, 100), (75, 100)}
    )

    mo.output.replace(
        mo.md(
            f"Consensus clustering:\n\n"
            f"| Variant | Clusters | Cophenetic Corr. |\n"
            f"|---------|----------|------------------|\n"
            f"| Full | {consensus_clst['cluster'].max()} | {coph_cor:.4f} |\n"
            f"| Filtered | {consensus_clst_filt['cluster'].max()} | {coph_cor_filt:.4f} |"
        )
    )
    return (
        df_consensus,
        consensus_clst,
        link,
        coph_cor,
        df_consensus_filt,
        consensus_clst_filt,
        link_filt,
        coph_cor_filt,
    )


@app.cell
def _():
    mo.md(
        """
        ## Consensus Matrix (Filtered)

        The filtered consensus excludes the (50–100 %) and (75–100 %)
        submatrices, which overlap heavily with core genes and can dilute
        the accessory-genome signal. Both full and filtered clustermaps are
        plotted below.
        """
    )


@app.cell
def _(FIG, consensus_clst, consensus_clst_filt, df_consensus, df_consensus_filt, link, link_filt):
    """Plot consensus clustermaps and cluster size bar charts."""
    # --- Clustermap: Full consensus ---
    g_full = sns.clustermap(
        df_consensus,
        figsize=(9, 9),
        row_linkage=link,
        col_linkage=link,
        yticklabels=False,
        xticklabels=False,
        cmap="hot_r",
    )
    g_full.fig.suptitle("Consensus Matrix (All Submatrices)", y=1.02)
    g_full.savefig(os.path.join(FIG, "4a_consensus_clustermap.png"), bbox_inches="tight")

    # --- Clustermap: Filtered consensus ---
    g_filt = sns.clustermap(
        df_consensus_filt,
        figsize=(9, 9),
        row_linkage=link_filt,
        col_linkage=link_filt,
        yticklabels=False,
        xticklabels=False,
        cmap="hot_r",
    )
    g_filt.fig.suptitle("Consensus Matrix (Filtered: excl. (50,100) & (75,100))", y=1.02)
    g_filt.savefig(os.path.join(FIG, "4a_consensus_clustermap_filtered.png"), bbox_inches="tight")

    # --- Bar charts: Cluster sizes ---
    fig_bars, (ax_bar1, ax_bar2) = plt.subplots(1, 2, figsize=(14, 4))

    counts_full = consensus_clst["cluster"].value_counts().sort_index()
    ax_bar1.bar(counts_full.index, counts_full.values)
    ax_bar1.set_xlabel("Cluster")
    ax_bar1.set_ylabel("Number of Genomes")
    ax_bar1.set_title(f"Full Consensus ({consensus_clst['cluster'].max()} clusters)")

    counts_filt = consensus_clst_filt["cluster"].value_counts().sort_index()
    ax_bar2.bar(counts_filt.index, counts_filt.values)
    ax_bar2.set_xlabel("Cluster")
    ax_bar2.set_ylabel("Number of Genomes")
    ax_bar2.set_title(f"Filtered Consensus ({consensus_clst_filt['cluster'].max()} clusters)")

    plt.tight_layout()
    fig_bars.savefig(os.path.join(FIG, "4a_consensus_cluster_sizes.png"), bbox_inches="tight")

    mo.output.replace(mo.vstack([g_full.fig, g_filt.fig, fig_bars]))


@app.cell
def _():
    mo.md(
        """
        ## NMF Output Shapes

        Final matrices saved for downstream analysis:

        - **L** (genomes × phylons): phylon membership weights per genome
        - **A** (phylons × genes): gene contribution weights per phylon
        - **L_binarized / A_binarized**: thresholded binary assignments
        - **Consensus clusters**: Ward-clustered genome groups from filtered consensus
        """
    )


@app.cell
def _(
    DATA,
    OUT,
    best_L_binarized_dict,
    best_L_norm_dict,
    best_A_binarized_dict,
    best_A_norm_dict,
    best_ranks_dict,
    coph_cor,
    coph_cor_filt,
    consensus_clst,
    consensus_clst_filt,
    df_metrics,
):
    """Save NMF outputs, consensus clusters, and summary CSV."""
    _best_rank = best_ranks_dict["full"]
    _nmf_dir = os.path.join(DATA, "processed", "nmf-outputs")

    # NMF matrices (index=True — genes/genomes are natural indices)
    best_L_norm_dict["full"].to_csv(os.path.join(_nmf_dir, "L.csv"))
    best_A_norm_dict["full"].to_csv(os.path.join(_nmf_dir, "A.csv"))
    best_L_binarized_dict["full"].to_csv(os.path.join(_nmf_dir, "L_binarized.csv"))
    best_A_binarized_dict["full"].to_csv(os.path.join(_nmf_dir, "A_binarized.csv"))

    # Consensus cluster assignments
    _consensus_out = consensus_clst_filt.reset_index()
    _consensus_out.columns = ["genome_id", "cluster"]
    _consensus_out.to_csv(os.path.join(OUT, "data", "4a_consensus_clusters.csv"), index=False)

    # Summary CSV
    _best_aic = df_metrics.loc[_best_rank, "AIC"]
    _summary = pd.DataFrame(
        {
            "metric": [
                "best_rank",
                "best_aic",
                "n_clusters_full",
                "cophenetic_full",
                "n_clusters_filtered",
                "cophenetic_filtered",
            ],
            "value": [
                _best_rank,
                _best_aic,
                consensus_clst["cluster"].max(),
                round(coph_cor, 4),
                consensus_clst_filt["cluster"].max(),
                round(coph_cor_filt, 4),
            ],
        }
    )
    _summary.to_csv(os.path.join(OUT, "data", "4a_nmf_summary.csv"), index=False)

    mo.output.replace(
        mo.vstack([
            mo.md(
                f"Saved NMF outputs:\n\n"
                f"- `{_nmf_dir}/L.csv` — L_norm ({best_L_norm_dict['full'].shape})\n"
                f"- `{_nmf_dir}/A.csv` — A_norm ({best_A_norm_dict['full'].shape})\n"
                f"- `{_nmf_dir}/L_binarized.csv` ({best_L_binarized_dict['full'].shape})\n"
                f"- `{_nmf_dir}/A_binarized.csv` ({best_A_binarized_dict['full'].shape})\n"
                f"- `{os.path.join(OUT, 'data')}/4a_consensus_clusters.csv` ({_consensus_out.shape[0]} genomes)\n"
                f"- `{os.path.join(OUT, 'data')}/4a_nmf_summary.csv`\n\n"
                f"Best rank: **{_best_rank}** (AIC: {_best_aic:.1f})"
            ),
            mo.ui.table(_summary),
        ])
    )


if __name__ == "__main__":
    app.run()
