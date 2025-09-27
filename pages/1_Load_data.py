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
# 2. Upload 10X files individually
elif option == "Upload 10X files (matrix + genes/features + barcodes)":
    matrix_file   = st.file_uploader("Upload matrix.mtx / matrix.mtx.gz", type=["mtx", "gz"])
    genes_file    = st.file_uploader("Upload genes.tsv / features.tsv.gz", type=["tsv", "gz"])
    barcodes_file = st.file_uploader("Upload barcodes.tsv / barcodes.tsv.gz", type=["tsv", "gz"])

    if matrix_file and genes_file and barcodes_file:
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                # Save uploaded files with expected names
                with open(os.path.join(tmpdir, "matrix.mtx.gz" if matrix_file.name.endswith(".gz") else "matrix.mtx"), "wb") as f:
                    f.write(matrix_file.read())
                with open(os.path.join(tmpdir, "features.tsv.gz" if genes_file.name.endswith(".gz") else "genes.tsv"), "wb") as f:
                    f.write(genes_file.read())
                with open(os.path.join(tmpdir, "barcodes.tsv.gz" if barcodes_file.name.endswith(".gz") else "barcodes.tsv"), "wb") as f:
                    f.write(barcodes_file.read())

                # Load with Scanpy
                adata = sc.read_10x_mtx(tmpdir, var_names="gene_symbols", cache=True)
                st.success("✅ Data loaded from 10X files!")

                # Save AnnData to temporary .h5ad file
                tmp_h5ad = os.path.join(tmpdir, "uploaded_data.h5ad")
                adata.write_h5ad(tmp_h5ad)

                # Read file back into memory for download
                with open(tmp_h5ad, "rb") as f:
                    h5ad_bytes = f.read()

                # Provide download button
                st.download_button(
                    label="💾 Download as .h5ad",
                    data=h5ad_bytes,
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

    # Summary
    st.subheader("🔍 Preview")
    st.write(f"**Number of cells:** {adata.n_obs}")
    st.write(f"**Number of genes:** {adata.n_vars}")

    try:
        df = pd.DataFrame(
            adata.X[:5, :5].toarray() if hasattr(adata.X, "toarray") else adata.X[:5, :5],
            index=adata.obs_names[:5],
            columns=adata.var_names[:5]
        )
        st.dataframe(df)
    except Exception as e:
        st.error(f"Could not display preview matrix: {e}")
else:
    st.info("Please upload a dataset or select demo data to continue.")
