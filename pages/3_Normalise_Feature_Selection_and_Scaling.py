import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile, os

st.header("Step 3: Normalisation, Feature Selection, and Scaling")

# --- Check if adata exists from Step 2 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please load and filter data first (Steps 1 & 2).")
    st.stop()

adata = st.session_state["adata"]

# --- Normalization ---
st.subheader("Normalization")

scale_factor = st.number_input(
    "Scale factor (target counts per cell)", 
    min_value=1000, value=10000, step=1000
)

if st.button("Run Normalisation", key="normalize"):
    sc.pp.normalize_total(adata, target_sum=scale_factor)
    sc.pp.log1p(adata)

    # Save log-normalized data as .raw for downstream DE / CellTypist
    adata.raw = adata.copy()

    st.success(f"✅ Normalized each cell to {scale_factor} counts and log-transformed.")
    st.session_state["adata"] = adata
    st.session_state["normalized"] = True

# --- Highly Variable Genes ---
st.subheader("Highly Variable Genes (HVGs)")

n_top_genes = st.number_input(
    "Number of variable genes", 
    min_value=500, value=2000, step=500
)

if st.button("Identify HVGs", key="hvg"):
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        flavor="seurat"
    )
    st.session_state["adata"] = adata
    st.session_state["hvg_done"] = True
    st.success(f"✅ Identified top {n_top_genes} highly variable genes.")

# --- Show HVG plots if done ---
if st.session_state.get("hvg_done", False):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: normalized dispersion
    axes[0].scatter(
        adata.var["means"][~adata.var["highly_variable"]],
        adata.var["dispersions_norm"][~adata.var["highly_variable"]],
        c="black", s=5, label="Other genes"
    )
    axes[0].scatter(
        adata.var["means"][adata.var["highly_variable"]],
        adata.var["dispersions_norm"][adata.var["highly_variable"]],
        c="red", s=5, label="Highly variable genes"
    )
    axes[0].set_xlabel("Mean expression of genes")
    axes[0].set_ylabel("Dispersions (normalized)")
    axes[0].set_title("Normalized")
    axes[0].legend(frameon=False)

    # Right: unnormalized dispersion
    axes[1].scatter(
        adata.var["means"][~adata.var["highly_variable"]],
        adata.var["dispersions"][~adata.var["highly_variable"]],
        c="black", s=5
    )
    axes[1].scatter(
        adata.var["means"][adata.var["highly_variable"]],
        adata.var["dispersions"][adata.var["highly_variable"]],
        c="red", s=5
    )
    axes[1].set_xlabel("Mean expression of genes")
    axes[1].set_ylabel("Dispersions (not normalized)")
    axes[1].set_title("Unnormalized")

    st.pyplot(fig)

# --- Scaling (all genes) ---
st.subheader("Scaling")

if st.button("Run Scaling", key="scale"):
    # Scale ALL genes (HVGs are still flagged for PCA/clustering)
    sc.pp.scale(adata, max_value=10)

    st.success("✅ Scaled all genes to unit variance and mean 0.")
               

    st.session_state["adata"] = adata

    # --- Download scaled AnnData ---
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        adata.write(tmp.name)
        tmp_path = tmp.name

    with open(tmp_path, "rb") as f:
        st.download_button(
            label="💾 Download normalized & scaled data (.h5ad)",
            data=f,
            file_name="adata_scaled.h5ad",
            mime="application/octet-stream",
            key="download_scaled"
        )
    os.remove(tmp_path)
