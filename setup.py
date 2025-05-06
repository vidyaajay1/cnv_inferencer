import os
from setuptools import setup, find_packages

setup(
    name="cnv_inferencer",  # or whatever your package is called
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
    "scanpy",
    "anndata",
    "pandas",
    "numpy",
    "mygene",
    "matplotlib",
    "scipy",
    "seaborn",
    "igraph",
    "leidenalg"
    ],

    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    python_requires=">=3.7",
)
