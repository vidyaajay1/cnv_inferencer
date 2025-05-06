import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse
from collections import defaultdict

MAIN_CHROMS = list(map(str, range(1, 23))) + ["X", "Y"]
WINDOW_GENES = 150
GAIN_THRES = +0.5
LOSS_THRES = -0.5
MIN_CELLS_SEG = 30

def call_cnvs(adata, output_dir="files", step=200, plot=True):
    """
    Calls CNVs from an AnnData object and saves summary files to disk.
    Returns the updated AnnData with CNV calls in .obs['cnv_calls'].
    """
    os.makedirs(output_dir, exist_ok=True)

    # 1. Filter for valid chromosomes
    mask = (
        adata.var["chromosome"].astype(str).isin(MAIN_CHROMS) &
        adata.var["start"].notna() & adata.var["end"].notna()
    )
    adata = adata[:, mask].copy()
    adata.var["chrom"] = adata.var["chromosome"].astype(str)
    adata.var["chrom"] = pd.Categorical(adata.var["chrom"], categories=MAIN_CHROMS, ordered=True)
    adata.var["start"] = adata.var["start"].astype(int)
    adata.var.sort_values(["chrom", "start"], inplace=True)
    adata = adata[:, adata.var.index].copy()

    # 2. Normalize and log2-transform (CPM)
    counts = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    libsize = counts.sum(1)
    normed = counts / libsize[:, None] * 1e6
    log2expr = np.log2(normed + 1)
    log2_df = pd.DataFrame(log2expr, index=adata.obs_names, columns=adata.var_names)

    # 3. Smoothing
    smooth_list = []
    for chrom in MAIN_CHROMS:
        g_chr = adata.var_names[adata.var["chrom"] == chrom]
        if len(g_chr) == 0:
            continue
        smooth_chr = log2_df[g_chr].rolling(window=WINDOW_GENES, axis=1, min_periods=1, center=True).mean()
        smooth_list.append(smooth_chr)
    smoothed_df = pd.concat(smooth_list, axis=1)

    # 4. Diploid reference cluster
    cnv_adata = sc.AnnData(X=smoothed_df.values, obs=adata.obs.copy(), var=pd.DataFrame(index=smoothed_df.columns))
    sc.pp.pca(cnv_adata, n_comps=30)
    sc.pp.neighbors(cnv_adata, n_neighbors=15, metric="correlation")
    sc.tl.leiden(cnv_adata, resolution=0.1, key_added="cnv_leiden")
    cnv_adata.obs["cnv_burden"] = np.mean(np.abs(cnv_adata.X), axis=1)
    ref_cluster = cnv_adata.obs.groupby("cnv_leiden")["cnv_burden"].median().idxmin()
    ref_cells = cnv_adata.obs_names[cnv_adata.obs["cnv_leiden"] == ref_cluster]
    ref_profile = smoothed_df.loc[ref_cells].median(axis=0)

    # 5. Relative CNV signal
    rel_df = smoothed_df.subtract(ref_profile, axis=1)

    # 6. Heatmap
    if plot:
        cluster_order = cnv_adata.obs.groupby("cnv_leiden")["cnv_burden"].median().sort_values().index
        ordered_cells, cluster_sizes = [], []
        for cl in cluster_order:
            members = cnv_adata.obs_names[cnv_adata.obs["cnv_leiden"] == cl]
            ordered_cells.extend(members)
            cluster_sizes.append(len(members))
        plot_df = rel_df.loc[ordered_cells]
        chroms = adata.var["chrom"].tolist()
        chrom_edges, chrom_labels = [0], []
        for i in range(1, len(chroms)):
            if chroms[i] != chroms[i-1]:
                chrom_edges.append(i)
                chrom_labels.append(chroms[i-1])
        chrom_edges.append(len(chroms))
        chrom_labels.append(chroms[-1])
        plot_mat = plot_df.iloc[:, ::step]

        fig, ax = plt.subplots(figsize=(16, 6))
        sns.heatmap(plot_mat, cmap="bwr", vmin=-0.6, vmax=0.6, xticklabels=False, yticklabels=False,
                    ax=ax, cbar_kws=dict(label="log2 FC vs. diploid ref"))
        cum_rows = np.cumsum(cluster_sizes)
        for y in cum_rows[:-1]: ax.axhline(y, color="k", lw=1.0)
        for x in chrom_edges: ax.axvline(x / step, color="k", lw=1.0)
        xtick_pos = [(chrom_edges[i] + chrom_edges[i+1]) / (2 * step) for i in range(len(chrom_labels))]
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(chrom_labels, rotation=90, fontsize=8)
        ytick_pos = [(cum_rows[i-1] if i else 0) + size/2 for i, size in enumerate(cluster_sizes)]
        ylabels = [f"{cl} (ref)" if cl == ref_cluster else str(cl) for cl in cluster_order]
        ax.set_yticks(ytick_pos)
        ax.set_yticklabels(ylabels, fontsize=8)
        ax.set_title("Smoothed CNV signal by chromosome and cnv_leiden")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "cnv_heatmap.png"))
        plt.close()

    # 7. Consensus segmentation
    cnv_cells = cnv_adata.obs_names[cnv_adata.obs["cnv_leiden"] != ref_cluster]
    cons_med = rel_df.loc[cnv_cells].median(axis=0)
    segments = []
    for chrom in MAIN_CHROMS:
        genes_idx = adata.var.index[adata.var["chrom"] == chrom]
        if len(genes_idx) == 0: continue
        sig = cons_med.loc[genes_idx].values
        pos = adata.var.loc[genes_idx, "start"].values
        cls = np.where(sig > GAIN_THRES, 1, np.where(sig < LOSS_THRES, -1, 0))
        run_start = 0
        for i in range(1, len(cls) + 1):
            if i == len(cls) or cls[i] != cls[run_start]:
                state = cls[run_start]
                if state != 0:
                    segments.append(dict(
                        chrom=chrom,
                        g0=genes_idx[run_start],
                        g1=genes_idx[i-1],
                        start=int(pos[run_start]),
                        end=int(pos[i-1]),
                        CN_call=2 + state
                    ))
                run_start = i

    # 8. Per-cell CNV call
    records = []
    for seg in segments:
        i0 = adata.var.index.get_loc(seg["g0"])
        i1 = adata.var.index.get_loc(seg["g1"]) + 1
        seg_vals = rel_df.iloc[:, slice(i0, i1)].mean(axis=1)
        cells_with = seg_vals > GAIN_THRES if seg["CN_call"] > 2 else seg_vals < LOSS_THRES
        if int(cells_with.sum()) >= MIN_CELLS_SEG:
            records.append(dict(
                segment=f"{seg['chrom']}:{seg['start']}-{seg['end']} (CN {seg['CN_call']})",
                n_cells=int(cells_with.sum())
            ))

    pd.DataFrame(records).sort_values("n_cells", ascending=False).to_csv(
        os.path.join(output_dir, "cnv_segments_summary.csv"), index=False)

    # Cell-level calls
    cells_with_cnv = []
    for seg in segments:
        i0 = adata.var.index.get_loc(seg["g0"])
        i1 = adata.var.index.get_loc(seg["g1"]) + 1
        seg_vals = rel_df.iloc[:, slice(i0, i1)].mean(axis=1)
        cells_with = seg_vals > GAIN_THRES if seg["CN_call"] > 2 else seg_vals < LOSS_THRES
        for cell, has_cnv in zip(adata.obs_names, cells_with):
            if has_cnv:
                cells_with_cnv.append({
                    "cell": cell,
                    "segment": f"{seg['chrom']}:{seg['start']}-{seg['end']} (CN {seg['CN_call']})",
                    "CNV_call": seg["CN_call"]
                })
    df_cells = pd.DataFrame(cells_with_cnv)
    df_cells.to_csv(os.path.join(output_dir, "cells_with_cnv.csv"), index=False)

    # Annotate AnnData with CNV calls
    cell2segs = defaultdict(list)
    for row in cells_with_cnv:
        cell2segs[row["cell"]].append(row["segment"])
    adata.obs["cnv_calls"] = (
        pd.Series({c: "; ".join(sorted(v)) for c, v in cell2segs.items()})
          .reindex(adata.obs_names).fillna("").astype("string")
    )
    adata.write_h5ad(os.path.join(output_dir, "sim_CNV_final_calls.h5ad"))
    print(f"CNV results saved to {output_dir}/")

    return adata
