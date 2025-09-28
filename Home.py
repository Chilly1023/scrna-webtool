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

You can explore your single-cell dataset in **six main steps**:

1. **📂 Load Data** – Upload `.h5ad` files, raw 10X input, or use the included PBMC3k demo dataset.  
2. **🔧 Preprocessing** – Perform quality control (QC), normalization, selection of highly variable genes (HVGs), and scaling.  
3. **📉 Linear Dimensional Reduction (PCA)** – Reduce dimensionality to highlight major sources of variation and prepare for clustering.  
4. **🔗 Clustering & UMAP** – Group cells into clusters (Leiden algorithm) and visualize them in 2D with UMAP.  
5. **🧬 Differential Expression (Marker Genes)** – Identify genes that distinguish clusters.  
6. **🧭 Gene Expression & Cell Type Annotation** – Visualize gene expression on UMAP, automatically detect cluster markers, and annotate clusters using CellTypist.  

---

## 📦 Data Input Options
- Upload `.h5ad` file (recommended for large datasets).  
- Upload **raw 10X files** (`matrix.mtx`, `genes.tsv/features.tsv`, `barcodes.tsv`) → auto-converts to `.h5ad`.  
- Use included **PBMC3k demo dataset**.  

---

ℹ️ This project was built with the help of **ChatGPT-5 (OpenAI)** for code structure, deployment, and documentation.
""")
