# CNV Inferencer

`cnv_inferencer` is a Python package for preprocessing single-cell RNA-seq data and inferring **copy number variations (CNVs)** at the cell level. It supports end-to-end workflows including quality control, clustering, annotation, CNV signal smoothing, reference cluster detection, and per-cell CNV calling with output visualizations and summary CSVs.

---

## Installation

Clone the repository and install in **editable mode**:

```bash
git clone https://github.com/vidyaajay1/cnv_inferencer.git
cd cnv_inferencer
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```
## Dependencies

Required Python packages:

scanpy \
anndata \
pandas \
numpy \
mygene \
matplotlib \
scipy \
seaborn \
igraph \
leidenalg 

## Usage

Sample Pipeline

```python
import scanpy as sc
from cnv_inferencer import (
    preprocess_and_cluster,
    annotate_genomic_positions,
    call_cnvs
)

# Load simulated or experimental scRNA-seq data
adata = sc.read_h5ad("files/PBMC_simCNV_2.h5ad")

# Preprocess & cluster
adata = preprocess_and_cluster(adata)

# Add genomic annotations (if missing)
adata = annotate_genomic_positions(adata)

# Call CNVs and save output
adata = call_cnvs(adata, output_dir="files")

# View results
print(adata.obs["cnv_calls"].head())
```

## Outputs

Calling CNVs with call_cnvs() will generate: 

files/cnv_heatmap.png — smoothed CNV heatmap \
files/cnv_segments_summary.csv — CNV segment counts \
files/cells_with_cnv.csv — per-cell CNV calls \
files/sim_CNV_final_calls.h5ad — updated AnnData object with `adata.obs["cnv_calls"]`

## Functions

preprocess_and_cluster(adata):	Performs QC, normalization, clustering, and UMAP \
annotate_genomic_positions(adata):	Queries gene coordinates via MyGene.info and adds chromosome, start, end to adata.var \
call_cnvs(adata, output_dir="files"):	Infers CNVs, saves heatmaps, CSVs, and updates `adata.obs["cnv_calls"]`

## Sample Data

If you’re testing with simulated CNVs, make sure `adata.var` contains gene IDs as gene_ids and that MyGene.info annotations can be fetched for those IDs. \
The package expects a `.h5ad` file that has cell_type in `adata.obs` and gene_ids in `adata.var`.\
Place your `.h5ad` test files here to run the CNV inference pipeline locally.

## What to look for in a good run of the package
Some warnings will be displayed on a good run, they're harmless and help in detecting the status of the code and debugging. For reference, a normal run
will look like this in the terminal window:
```
(.venv) vidyaajay@Vidyas-MacBook-Pro cnv_inferencer % python test_run.py
Input sequence provided is already in string format. No operation performed
Input sequence provided is already in string format. No operation performed
5 input query terms found dup hits:	[('ENSG00000234162', 2), ('ENSG00000227110', 2), ('ENSG00000249738', 2), ('ENSG00000280018', 2), ('E
411 input query terms found no hit:	['ENSG00000238009', 'ENSG00000230699', 'ENSG00000236948', 'ENSG00000277726', 'ENSG00000271895', 'ENS
Input sequence provided is already in string format. No operation performed
Input sequence provided is already in string format. No operation performed
5 input query terms found dup hits:	[('ENSG00000234162', 2), ('ENSG00000227110', 2), ('ENSG00000249738', 2), ('ENSG00000280018', 2), ('E
395 input query terms found no hit:	['ENSG00000238009', 'ENSG00000230699', 'ENSG00000236948', 'ENSG00000277726', 'ENSG00000271895', 'ENS
Index(['gene_ids', 'feature_types', 'genome', 'mt', 'ribo',
       'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts',
       'total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions',
       'dispersions_norm', 'chromosome', 'start', 'end', 'strand'],
      dtype='object')
/Users/vidyaajay/cnv_inferencer/cnv_inferencer/caller.py:57: FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.

 To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.  directed must also be False to work with igraph's implementation.
  sc.tl.leiden(cnv_adata, resolution=0.1, key_added="cnv_leiden")
/Users/vidyaajay/cnv_inferencer/cnv_inferencer/caller.py:59: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
  ref_cluster = cnv_adata.obs.groupby("cnv_leiden")["cnv_burden"].median().idxmin()
/Users/vidyaajay/cnv_inferencer/cnv_inferencer/caller.py:68: FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas. Pass observed=False to retain current behavior or observed=True to adopt the future default and silence this warning.
  cluster_order = cnv_adata.obs.groupby("cnv_leiden")["cnv_burden"].median().sort_values().index
CNV results saved to files/
AAACCCAAGCGCCCAT-1                                                     
AAACCCAAGGTTCCGC-1    10:101845599-102065349 (CN 3); 10:103453240-10...
AAACCCACAGAGTTGG-1    12:7089587-7115736 (CN 3); 19:57435325-5746666...
AAACCCACAGGTATGG-1    10:103453240-104309698 (CN 3); 10:110868890-11...
AAACCCACATAGTCAC-1                                                     
Name: cnv_calls, dtype: category
Categories (6444, object): ['', '1:10456522-10456522 (CN 3)',
                            '1:10456522-10456522 (CN 3); 1:10472288-173081..., '1:10456522-10456522 (CN 3); 1:10472288-173081...,
                            ..., 'X:54192823-54530183 (CN 3)',
                            'X:72301638-73943934 (CN 3)', 'X:72301638-73943934 (CN 3); X:76173040-775048...,
                            'X:76173040-77504880 (CN 3)']
```

## Questions

For any questions or help in running the package please contact vajay1@jh.edu and I'll get back to you quickly!
