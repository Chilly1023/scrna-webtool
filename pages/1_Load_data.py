import streamlit as st
import scanpy as sc
import tempfile, os, io
import pandas as pd

# st.set_page_config(page_title="Load Data", page_icon="📂")

# --- Sidebar controls ---
st.sidebar.header("Step 1: Load Data 📂")

option = st.sidebar.radio(
    "Choose how to load your data:",
    ["Upload files (.h5ad)", 
     "Upload 10X files (matrix + genes/features + barcodes)", 
     "Use Demo Data"],
    index=2
)

adata = None

# 1. Upload h5ad
if option == "Upload files (.h5ad)":
    uploaded_files = st.sidebar.file_uploader(
        "Upload one or more .h5ad files",
        type=["h5ad"],
        accept_multiple_files=True
    )
    if uploaded_files:
        adata_list = []
        for f in uploaded_files:
            try:
                adata_temp = sc.read_h5ad(io.BytesIO(f.read()))
                adata_list.append(adata_temp)
                st.sidebar.success(f"✅ Loaded {f.name} successfully!")
            except Exception as e:
                st.sidebar.error(f"❌ Failed to load {f.name}: {e}")

        # Merge if multiple
        try:
            if len(adata_list) == 1:
                adata = adata_list[0]
            else:
                adata = adata_list[0].concatenate(*adata_list[1:], batch_key="batch")
                st.sidebar.success("✅ Multiple .h5ad files merged into one AnnData object with `batch` column.")
        except Exception as e:
            st.sidebar.error(f"❌ Failed to merge .h5ad files: {e}")


# 2. Upload 10X files individually
elif option == "Upload 10X files (matrix + genes/features + barcodes)":
    uploaded_files = st.sidebar.file_uploader(
        "Upload your 10X files (matrix.mtx + genes.tsv/features.tsv + barcodes.tsv)",
        type=["mtx", "tsv", "gz"],
        accept_multiple_files=True
    )

    if uploaded_files:
        # classify by name
        mtx_file = next((f for f in uploaded_files if f.name.endswith(".mtx") or f.name.endswith(".mtx.gz")), None)
        genes_file = next((f for f in uploaded_files if "genes" in f.name or "features" in f.name), None)
        barcodes_file = next((f for f in uploaded_files if "barcodes" in f.name), None)

        if mtx_file and genes_file and barcodes_file:
            with tempfile.TemporaryDirectory() as tmpdir:
                try:
                    # Decide target filenames
                    mtx_name = "matrix.mtx.gz" if mtx_file.name.endswith(".gz") else "matrix.mtx"
                    if genes_file.name.endswith(".gz"):
                        feat_name = "features.tsv.gz" if "features" in genes_file.name else "genes.tsv.gz"
                    else:
                        feat_name = "features.tsv" if "features" in genes_file.name else "genes.tsv"
                    bar_name = "barcodes.tsv.gz" if barcodes_file.name.endswith(".gz") else "barcodes.tsv"

                    # Save files
                    for f, fname in [(mtx_file, mtx_name), (genes_file, feat_name), (barcodes_file, bar_name)]:
                        with open(os.path.join(tmpdir, fname), "wb") as out_f:
                            out_f.write(f.read())

                    # Read with Scanpy
                    try:
                        adata = sc.read_10x_mtx(tmpdir, var_names="gene_symbols", cache=False)
                    except Exception:
                        adata = sc.read_10x_mtx(tmpdir, var_names="gene_ids", cache=False)

                    st.sidebar.success("✅ Data loaded from 10X files!")

                except Exception as e:
                    st.sidebar.error(f"❌ Error reading 10X files: {e}")

# 3. Use demo data (make relative path robust)
elif option == "Use Demo Data":
    try:
        adata = sc.datasets.pbmc3k()  # raw counts
        st.session_state["adata"] = adata
        st.sidebar.success("✅ Demo data loaded (PBMC 3k)!")
    except Exception as e:
        st.sidebar.error(f"❌ Could not load demo dataset: {e}")


# --- Main page explanation ---
st.title("📂 Load your data here!")

st.markdown("""
Welcome to the **data loading step** of the scRNA-seq Webtool!  
Here you can choose one of three options to bring your dataset into the analysis pipeline:

### 🔹 Option 1: Upload `.h5ad` files
- `.h5ad` is the native **AnnData format** used in Scanpy.  
- It can store not only the raw count matrix, but also preprocessing results (QC, normalization), embeddings (PCA/UMAP), clustering, and annotations.  
- If you upload a **single `.h5ad` file**, it will be loaded as is.  
- If you upload **multiple `.h5ad` files**, the tool will automatically **merge them** into one AnnData object and add a new column `batch` in `adata.obs` to track the origin of each cell.

### 🔹 Option 2: Upload 10X Genomics raw files
- Standard **Cell Ranger** outputs:  
  - `matrix.mtx` : sparse count matrix (genes × cells)  
  - `genes.tsv` or `features.tsv` : gene information  
  - `barcodes.tsv` : cell identifiers  
- Please select **all three files at once** (hold `Ctrl`/`Cmd` or `Shift` when selecting).  

### 🔹 Option 3: Use Demo Dataset
- Built-in **PBMC 3k dataset** (~2,700 peripheral blood mononuclear cells from 10X Genomics).  
- Perfect for learning and testing the workflow.  
""")

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
