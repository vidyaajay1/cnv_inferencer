import scanpy as sc
from cnv_inferencer import preprocess_and_cluster, annotate_genomic_positions, call_cnvs

adata = sc.read_h5ad("files/PBMC_simCNV_2.h5ad")
print(adata.var.columns)
adata = preprocess_and_cluster(adata)
adata = annotate_genomic_positions(adata)
print(adata.var.columns)
adata = call_cnvs(adata, output_dir="files")
print(adata.obs["cnv_calls"].head())
