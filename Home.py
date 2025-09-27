import streamlit as st

st.set_page_config(
    page_title="scRNA-seq Webtool",
    page_icon="🧬",
    layout="wide"
)

# Title
st.title("🧬 Single-cell RNA-seq Analysis Webtool")
st.markdown("---")

# Introduction
st.markdown("""
Welcome! 👋  

This webtool lets you process **single-cell RNA-seq data** step by step, 
similar to workflows in [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) (Python) 
and [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) (R).  

You can explore datasets interactively — from **quality control** to **clustering**, **dimensionality reduction**, 
and **marker gene discovery**.
""")

# Data loading options
st.subheader("📂 How to load your data")
st.markdown("""
You have three options to begin analysis:

1. **Upload a preprocessed file (.h5ad)**  
   - Fastest and most efficient.  
   - Recommended for large datasets.  

2. **Upload raw 10X files**  
   - Provide the three files individually:  
     - `matrix.mtx` (or `matrix.mtx.gz`)  
     - `genes.tsv` / `features.tsv.gz`  
     - `barcodes.tsv` / `barcodes.tsv.gz`  
   - These will be converted internally to `.h5ad`.  

3. **Use Demo Data**  
   - Loads a small PBMC3k dataset provided with the app (from your `data/` folder).  
   - Great for testing and learning.  
""")

# Workflow
st.subheader("🧪 Workflow")
st.markdown("""
Use the **sidebar navigation** to process your data step by step:

1. **Load Data**  
2. **QC Filtering**  
3. **Preprocessing** (normalization & highly variable genes)  
4. **Run PCA**  
5. **Clustering**  
6. **Run UMAP**  
7. **DEGs** (marker gene discovery)  
8. **Assign Cell Type Identity**  

---
""")

st.success("👉 Start by going to **Step 1: Load Data** in the sidebar.")
