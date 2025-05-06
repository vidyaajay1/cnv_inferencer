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

Calling CNVs with call_cnvs() will generate: \

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
