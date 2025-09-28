import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile, os
import pandas as pd
import numpy as np

st.title("🔧 Preprocessing")

st.markdown("""
This step prepares your single-cell data for downstream analysis.  
The preprocessing workflow includes:

1. **QC and Filtering** – assess sequencing quality and select cells for further analysis.  
2. **Normalization** – standardize sequencing depth across cells.  
3. **Highly Variable Gene (HVG) Selection** – identify informative features for dimensionality reduction and clustering.  
4. **Scaling** – shift gene expression to mean = 0 and variance = 1 so that all features contribute comparably to downstream dimensionality reduction methods.  
   
   *👉 It is recommended to perform these steps in the given order, as each step depends on the results of the previous one.*  
""")


# --- Check if adata exists from Step 1 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please load data in Step 1 first.")
    st.stop()

adata = st.session_state["adata"]

# --- Sidebar to choose submodule ---
submodule = st.sidebar.radio(
    "Choose preprocessing step:",
    ["QC and Filtering", "Normalization", "HVG Selection", "Scaling"]
)

def save_and_download(adata, filename, label="Download data (.h5ad)"):
    """Helper to save AnnData and provide a download button in Streamlit."""
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        adata.write(tmp.name)
        tmp_path = tmp.name
    with open(tmp_path, "rb") as f:
        st.download_button(
            label=f"💾 {label}",
            data=f,
            file_name=filename,
            mime="application/octet-stream"
        )
    os.remove(tmp_path)


# =========================================================
# --- QC and Filtering ---
# =========================================================
if submodule == "QC and Filtering":
    st.header("🧹 QC and Filtering")
    st.markdown("""
    In this step, you can evaluate cell quality and remove low-quality cells.  

    **Typical workflow:**
    1. Visualize QC metrics (before filtering).  
    2. Decide filtering thresholds (number of genes, total counts, mitochondrial %).  
    3. Apply filtering.  
    4. Re-plot QC metrics after filtering to check the effect.  
    """)

    # --- Initialize raw copy if not exists ---
    if "adata_raw" not in st.session_state:
        st.session_state["adata_raw"] = adata.copy()

    adata_raw = st.session_state["adata_raw"]

    # --- QC METRICS (before filtering) ---
    if "n_genes_by_counts" not in adata_raw.obs.columns:
        adata_raw.var["mt"] = adata_raw.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    st.subheader("QC Violin Plots (before filtering)")

    qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    qc_labels = {
        "n_genes_by_counts": "Number of genes",
        "total_counts": "Total counts",
        "pct_counts_mt": "Mitochondrial %"
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, metric in enumerate(qc_metrics):
        sc.pl.violin(
            adata_raw,
            keys=metric,
            groupby=None,
            jitter=0.4,
            multi_panel=False,
            ax=axes[i],
            show=False
        )
        axes[i].set_title(qc_labels[metric])
    st.pyplot(fig)

    # --- Scatter plots (user selects axes) ---
    # --- QC Scatter Plots ---
    st.subheader("QC Scatter Plots (before filtering)")

    # 反向映射（从友好标签 → 真实列名）
    qc_reverse = {v: k for k, v in qc_labels.items()}


    # Select axes for first scatter
    col1, col2 = st.columns(2)
    with col1:
        x1 = st.selectbox("Select X-axis (Plot 1)", list(qc_labels.values()), index=0)
        y1 = st.selectbox("Select Y-axis (Plot 1)", list(qc_labels.values()), index=2)

    with col2:
        x2 = st.selectbox("Select X-axis (Plot 2)", list(qc_labels.values()), index=1)
        y2 = st.selectbox("Select Y-axis (Plot 2)", list(qc_labels.values()), index=0)

    # Create figure with 2 subplots side by side
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1 (light red)
    axes[0].scatter(
        adata_raw.obs[qc_reverse[x1]],
        adata_raw.obs[qc_reverse[y1]],
        s=8,
        c="lightcoral",   # light red
        alpha=0.6,
        edgecolors="none"
    )
    axes[0].set_xlabel(x1)
    axes[0].set_ylabel(y1)
    axes[0].set_title(f"{y1} vs {x1}")

    # Plot 2 (light red)
    axes[1].scatter(
        adata_raw.obs[qc_reverse[x2]],
        adata_raw.obs[qc_reverse[y2]],
        s=8,
        c="lightcoral",   # light red
        alpha=0.6,
        edgecolors="none"
    )
    axes[1].set_xlabel(x2)
    axes[1].set_ylabel(y2)
    axes[1].set_title(f"{y2} vs {x2}")

    st.pyplot(fig)

    # --- Filtering thresholds ---
    st.subheader("Filtering thresholds")

    min_features = st.number_input("Minimum number of genes", min_value=0, value=200)
    max_features = st.number_input("Maximum number of genes", min_value=0, value=2500)
    max_percent_mt = st.number_input("Maximum mitochondrial percentage", min_value=0, value=5)

    if st.button("Apply filtering"):
        initial_cells = adata_raw.n_obs

        adata = adata_raw[
            (adata_raw.obs["n_genes_by_counts"] > min_features) &
            (adata_raw.obs["n_genes_by_counts"] < max_features) &
            (adata_raw.obs["pct_counts_mt"] < max_percent_mt),
            :
        ].copy()

        st.session_state["adata"] = adata
        st.success(f"Filtered from {initial_cells} cells to {adata.n_obs} cells.")

        # Save filtered AnnData
        save_and_download(adata, "filtered_data.h5ad", "Download filtered data (.h5ad)")

        # with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        #     adata.write(tmp.name)
        #     tmp_path = tmp.name
        # with open(tmp_path, "rb") as f:
        #     st.download_button("💾 Download filtered data (.h5ad)", f, "filtered_data.h5ad", "application/octet-stream")
        # os.remove(tmp_path)

        # --- QC plots after filtering ---
        st.subheader("QC Violin Plots (after filtering)")
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        for i, metric in enumerate(qc_metrics):
            sc.pl.violin(
                adata,
                keys=metric,
                groupby=None,
                jitter=0.4,
                multi_panel=False,
                ax=axes[i],
                show=False
            )
            axes[i].set_title(qc_labels[metric])
        st.pyplot(fig)


# =========================================================
# --- Normalization ---
# =========================================================
elif submodule == "Normalization":
    st.header("📊 Normalization")
    st.markdown("""
    Normalize sequencing depth across cells to make them comparable.  
    This ensures that differences reflect biology rather than sequencing depth.  
    """)

    scale_factor = st.number_input("Scale factor (target counts per cell)", min_value=1000, value=10000, step=1000)

    if st.button("Run Normalisation"):
        sc.pp.normalize_total(adata, target_sum=scale_factor)
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        st.success(f"✅ Normalized each cell to {scale_factor} counts and log-transformed.")
        st.session_state["adata"] = adata

        save_and_download(adata, "normalized_data.h5ad", "Download normalized data (.h5ad)")


# =========================================================
# --- HVG Selection ---
# =========================================================
elif submodule == "HVG Selection":
    st.header("✨ Highly Variable Gene (HVG) Selection")
    st.markdown("""
    Identify the most informative features that drive biological variability.  
    These HVGs are used for dimensionality reduction and clustering.  
    """)

    n_top_genes = st.number_input("Number of variable genes", min_value=500, value=2000, step=500)

    if st.button("Identify HVGs"):
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top_genes,
            min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat"
        )
        st.session_state["adata"] = adata
        st.session_state["hvg_done"] = True
        st.success(f"✅ Identified top {n_top_genes} highly variable genes.")
        save_and_download(adata, "hvg_data.h5ad", "Download HVG data (.h5ad)")


    if st.session_state.get("hvg_done", False):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
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
        axes[0].set_title("Normalized dispersion"); axes[0].legend(frameon=False)

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
        axes[1].set_title("Raw dispersion")
        st.pyplot(fig)

    

# =========================================================
# --- Scaling ---
# =========================================================
elif submodule == "Scaling":
    st.header("⚖️ Scaling")
    st.markdown("""
    Scale each feature to mean = 0 and variance = 1,  
    so that all features contribute comparably to downstream dimensionality reduction methods.  
    """)

    if st.button("Run Scaling"):
        sc.pp.scale(adata, max_value=10)
        st.success("✅ Scaled all genes to unit variance and mean 0.")
        st.session_state["adata"] = adata

        
        save_and_download(adata, "scaled_data.h5ad", "Download scaled data (.h5ad)")
