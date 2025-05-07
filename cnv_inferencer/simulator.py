import scanpy as sc
import numpy as np
import pandas as pd
import random
import scipy.sparse as sp

from .preprocessing import annotate_genomic_positions


def generate_cnv_region(chromosomes, min_size, max_size):
    chrom = random.choice(chromosomes)
    start = random.randint(1, 100_000_000)
    end = start + random.randint(min_size, max_size)
    return {
        'chromosome': chrom,
        'start': start,
        'end': end,
        'type': random.choice(['gain', 'loss'])
    }


def assign_cnvs_unique(adata, cell_types, num_cnvs, chromosomes, min_size, max_size, cnv_prob=0.1):
    cnvs = [generate_cnv_region(chromosomes, min_size, max_size) for _ in range(num_cnvs)]
    used = []
    assignments = {}

    for cell_type in cell_types:
        cells = adata.obs.query("cell_type == @cell_type").index.tolist()
        n_cnv = int(round(len(cells) * cnv_prob))
        if n_cnv == 0:
            assignments[cell_type] = {}
            continue

        selected = random.sample(cells, n_cnv)
        available = [c for c in cnvs if c not in used]
        if not available:
            break
        cnv = random.choice(available)
        used.append(cnv)

        assignments[cell_type] = {cell: {'cnv': cnv, 'type': cnv['type']} for cell in selected}

    return assignments


def simulate_cnv_impact(adata, cnv_assignments):
    counts = adata.layers['counts']
    for ct, mapping in cnv_assignments.items():
        for cell, info in mapping.items():
            cnv = info['cnv']
            mask = (
                (adata.var['chromosome'] == cnv['chromosome']) &
                (adata.var['start'] >= cnv['start']) &
                (adata.var['end'] <= cnv['end'])
            )
            genes = adata.var.index[mask]
            idxs = [adata.var_names.get_loc(g) for g in genes]

            factor = 3 if info['type'] == 'gain' else 0.5
            row_idx = adata.obs_names.get_loc(cell)

            for gi in idxs:
                counts[row_idx, gi] *= factor

    adata.layers['counts'] = counts
    return adata


def simulate_cnvs_on_adata(
    adata,
    num_cnvs=1,
    cnv_prob=0.1,
    min_size=100_000,
    max_size=500_000,
    seed=5
):
    """
    Adds simulated CNVs to the input AnnData object and returns a modified copy.
    """
    random.seed(seed)
    np.random.seed(seed)

    counts = adata.X
    if sp.issparse(counts):
        counts = counts.toarray()
    adata.layers['counts'] = counts.copy()

    adata = annotate_genomic_positions(adata)

    adata_cnv = adata.copy()
    cell_types = adata_cnv.obs['cell_type'].unique()
    chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

    assignments = assign_cnvs_unique(
        adata_cnv, cell_types, num_cnvs,
        chromosomes, min_size, max_size, cnv_prob
    )
    adata_cnv = simulate_cnv_impact(adata_cnv, assignments)

    adata_cnv.obs['simulated_cnvs'] = ''
    for ct, mapping in assignments.items():
        for cell, info in mapping.items():
            cnv = info['cnv']
            desc = f"{cnv['chromosome']}:{cnv['start']}-{cnv['end']} ({cnv['type']})"
            adata_cnv.obs.at[cell, 'simulated_cnvs'] = desc

    return adata_cnv
