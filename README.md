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