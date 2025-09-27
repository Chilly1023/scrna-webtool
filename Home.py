import streamlit as st

st.title("🧬 scRNA-seq Analysis Webtool")

st.markdown("""
Welcome to the **Single-cell RNA-seq Analysis Webtool** 👋  
This interactive app was developed as part of my **PhD coursework** to make
the [Scanpy](https://scanpy.readthedocs.io/) and [Seurat](https://satijalab.org/seurat/) 
workflows more accessible through a simple, step-by-step web interface.

The tool is inspired by:
- [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial) (Satija Lab)  
- [Scanpy PBMC3k tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)  

---

## 🚀 Workflow Overview

You can explore your single-cell dataset using the following steps:

1. **Load Data** – Upload `.h5ad` or raw 10X files, or use demo PBMC3k data  
2. **QC Filtering** – Remove low-quality cells and high-mitochondrial content  
3. **Normalisation & Feature Selection** – Log-normalisation and HVG selection  
4. **Linear Dimensional Reduction** – PCA to capture major variance  
5. **Clustering** – Leiden/Louvain clustering on the neighborhood graph  
6. **Non-linear Dimensional Reduction** – UMAP for 2D visualization  
7. **DEGs** – Identify marker genes for each cluster  
8. **Assign Cell Type Identity** – Annotate clusters based on known markers  

---

## 📦 Data Input Options
- Upload `.h5ad` file (recommended for large datasets)  
- Upload **raw 10X files** (`matrix.mtx`, `genes.tsv/features.tsv`, `barcodes.tsv`) → auto-converts to `.h5ad`  
- Use included **PBMC3k demo dataset**  

---

ℹ️ This project was built with the help of **ChatGPT-5 (OpenAI)** for code structure, deployment, and documentation.
""")
