import scanpy as sc
from cnv_inferencer import preprocess_and_cluster, annotate_genomic_positions, call_cnvs
from cnv_inferencer import simulate_cnvs_on_adata


adata = sc.read_h5ad("files/input_file.h5ad")
#adata_sim = simulate_cnvs_on_adata(adata, num_cnvs=3, cnv_prob=0.2)
#adata_sim.write("files/simulated_output.h5ad")
print(adata.var.columns)
adata = preprocess_and_cluster(adata)
adata = annotate_genomic_positions(adata)
print(adata.var.columns)
adata = call_cnvs(adata, output_dir="files")
print(adata.obs["cnv_calls"].head())
