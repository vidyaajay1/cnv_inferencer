import scanpy as sc
import numpy as np
import pandas as pd
import ast
import mygene
import scipy.sparse as sp
import matplotlib.pyplot as plt

def annotate_genomic_positions(adata):
    """
    Query MyGene.info to fetch genomic positions for each gene and merge into adata.var.
    """
    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        adata.var['gene_ids'].tolist(),
        scopes=['ensembl.gene', 'symbol'],
        fields='genomic_pos',
        species='human'
    )
    df = pd.DataFrame(results)

    def parse_pos(x):
        if isinstance(x, str):
            try:
                return ast.literal_eval(x)
            except Exception:
                return None
        if isinstance(x, list) and x:
            return x[0]
        if isinstance(x, dict):
            return x
        return None

    df['pos'] = df['genomic_pos'].apply(parse_pos)
    df = df.dropna(subset=['pos'])

    df[['chromosome', 'start', 'end', 'strand']] = df['pos'].apply(
        lambda d: pd.Series({
            'chromosome': d.get('chr'),
            'start': d.get('start'),
            'end': d.get('end'),
            'strand': d.get('strand')
        })
    )
    df = df.drop_duplicates(subset=['query'])

    # Remove any existing genomic columns before merging to prevent conflicts
    for col in ['chromosome', 'start', 'end', 'strand']:
        if col in adata.var.columns:
            adata.var.drop(columns=[col], inplace=True)

    # Set index for merge
    genomic_df = df[['query', 'chromosome', 'start', 'end', 'strand']].drop_duplicates()
    genomic_df.set_index('query', inplace=True)

    # Merge safely
    adata.var = adata.var.merge(
        genomic_df,
        left_on='gene_ids',
        right_index=True,
        how='left'
    )

    # Ensure proper types
    adata.var['chromosome'] = adata.var['chromosome'].astype(str)
    adata.var['start'] = pd.to_numeric(adata.var['start'], errors='coerce').astype('Int64')
    adata.var['end'] = pd.to_numeric(adata.var['end'], errors='coerce').astype('Int64')

    return adata


def preprocess_and_cluster(ad, resolution=0.5, n_neighbors=20, n_pcs=15):
    """
    Complete preprocessing and clustering pipeline:
      1. Add genomic annotation using MyGene.info
      2. QC metrics and filtering
      3. Normalize and log-transform
      4. Identify HVGs
      5. Dimensionality reduction and clustering

    Returns the annotated and clustered AnnData object.
    """
    ad = annotate_genomic_positions(ad)

    # Identify mitochondrial and ribosomal genes
    ad.var['mt'] = ad.var_names.astype(str).str.startswith('MT-')
    ribo_prefix = ("RPS", "RPL")
    ad.var['ribo'] = ad.var_names.astype(str).str.startswith(ribo_prefix)
    sc.pp.calculate_qc_metrics(ad, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

    # QC filtering
    ad = ad[ad.obs['pct_counts_mt'] < 20, :].copy()
    sc.pp.filter_cells(ad, min_genes=500)
    sc.pp.filter_cells(ad, max_counts=30000)
    sc.pp.filter_genes(ad, min_cells=3)

    # Normalize and log
    ad.layers['counts'] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.layers['lognorm'] = ad.X.copy()
    sc.pp.highly_variable_genes(ad, min_mean=0.0125, max_mean=6, min_disp=0.25)

    # PCA + clustering
    sc.tl.pca(ad, mask_var='highly_variable')
    sc.pp.neighbors(ad, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(ad, resolution=resolution)
    sc.tl.umap(ad)

    # Plot UMAP
    color_col = "cell_type" if "cell_type" in ad.obs.columns else "leiden"
    sc.pl.umap(ad, color=color_col, legend_loc="on data", frameon=False, show=False)
    plt.savefig('files/umap_clusters.png')
    plt.close()

    return ad
