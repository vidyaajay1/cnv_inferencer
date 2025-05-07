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

    # Ensure counts layer
    counts = adata.X
    if sp.issparse(counts):
        counts = counts.toarray()
    adata.layers['counts'] = counts.copy()

    # Add genomic positions
    adata = annotate_genomic_positions(adata)

    # Simulate CNVs
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
