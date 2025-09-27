import streamlit as st
import scanpy as sc
import tempfile, os, io
import pandas as pd

st.header("Step 1: Load Data")

option = st.radio(
    "Choose how to load your data:",
    ["Upload file (.h5ad)", 
     "Upload 10X files (matrix + genes/features + barcodes) and generate .h5ad file", 
     "Use Demo Data"],
    index=2
)

adata = None

# 1. Upload h5ad
if option == "Upload file (.h5ad)":
    uploaded_file = st.file_uploader("Upload your .h5ad file", type=["h5ad"])
    if uploaded_file:
        try:
            adata = sc.read_h5ad(uploaded_file)
            st.success("✅ Data loaded from .h5ad!")
        except Exception as e:
            st.error(f"❌ Error reading .h5ad: {e}")

# 2. Upload 10X files individually
elif option == "Upload 10X files (matrix + genes/features + barcodes) and generate .h5ad file":
    matrix_file   = st.file_uploader("Upload matrix.mtx / matrix.mtx.gz", type=["mtx", "gz"])
    genes_file    = st.file_uploader("Upload genes.tsv / features.tsv(.gz)", type=["tsv", "gz"])
    barcodes_file = st.file_uploader("Upload barcodes.tsv / barcodes.tsv.gz", type=["tsv", "gz"])

    if matrix_file and genes_file and barcodes_file:
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                # Decide target filenames  (support v2 genes.tsv or v3 features.tsv)
                mtx_name = "matrix.mtx.gz" if matrix_file.name.endswith(".gz") else "matrix.mtx"
                # If user provided features.tsv(.gz) keep that name; if genes.tsv keep genes.tsv
                if genes_file.name.endswith(".gz"):
                    feat_name = "features.tsv.gz" if "features" in genes_file.name else "genes.tsv.gz"
                else:
                    feat_name = "features.tsv" if "features" in genes_file.name else "genes.tsv"
                bar_name = "barcodes.tsv.gz" if barcodes_file.name.endswith(".gz") else "barcodes.tsv"

                # Save uploaded files into temp dir with expected names
                with open(os.path.join(tmpdir, mtx_name), "wb") as f:
                    f.write(matrix_file.read())
                with open(os.path.join(tmpdir, feat_name), "wb") as f:
                    f.write(genes_file.read())
                with open(os.path.join(tmpdir, bar_name), "wb") as f:
                    f.write(barcodes_file.read())

                # Try reading with gene symbols first; if it fails, fall back to gene_ids
                try:
                    adata = sc.read_10x_mtx(tmpdir, var_names="gene_symbols", cache=False)
                except Exception:
                    adata = sc.read_10x_mtx(tmpdir, var_names="gene_ids", cache=False)

                st.success("✅ Data loaded from 10X files!")

                # Save AnnData object into memory buffer for download
                buffer = io.BytesIO()
                adata.write_h5ad(buffer)
                buffer.seek(0)

                st.download_button(
                    label="💾 Download as .h5ad",
                    data=buffer,
                    file_name="uploaded_data.h5ad",
                    mime="application/octet-stream"
                )

            except Exception as e:
                st.error(f"❌ Error reading 10X files: {e}")

# 3. Use demo data (your raw 10X files in data/)
elif option == "Use Demo Data":
    try:
        adata = sc.read_10x_mtx("data", var_names="gene_symbols", cache=True)
        st.success("✅ Demo data loaded!")
    except Exception as e:
        st.error(f"❌ Could not load demo dataset: {e}")

# Store and preview
if adata is not None:
    st.session_state.adata = adata

    st.subheader("🔍 Preview")
    st.write(f"**Number of cells:** {adata.n_obs}")
    st.write(f"**Number of genes:** {adata.n_vars}")

    try:
        import numpy as np
        X = adata.X[:5, :5].toarray() if hasattr(adata.X, "toarray") else adata.X[:5, :5]
        df = pd.DataFrame(X, index=adata.obs_names[:5], columns=adata.var_names[:5])
        st.dataframe(df)
    except Exception as e:
        st.error(f"Could not display preview matrix: {e}")
else:
    st.info("Please upload a dataset or select demo data to continue.")
