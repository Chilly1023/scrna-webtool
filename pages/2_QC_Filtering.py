import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import tempfile
import os

st.header("Step 2: QC and Filtering")

# --- Check if adata exists from Step 1 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please load data in Step 1 first.")
    st.stop()

adata = st.session_state["adata"]

# --- Calculate QC metrics (same as Scanpy tutorial) ---
if "n_genes_by_counts" not in adata.obs.columns:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")  # annotate mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# --- Violin plots (genes, counts, mt%) ---
st.subheader("QC Violin Plots")

qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
qc_labels = {
    "n_genes_by_counts": "Number of genes",
    "total_counts": "Total counts",
    "pct_counts_mt": "Mitochondrial %"
}

fig, axes = plt.subplots(1, 3, figsize=(15, 4))  # 1 row, 3 columns

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

# --- Scatter plots (user selects axes) ---
# --- QC Scatter Plots ---
st.subheader("QC Scatter Plots")

qc_metrics = {
    "Number of genes": "n_genes_by_counts",
    "Total counts": "total_counts",
    "Mitochondrial %": "pct_counts_mt"
}

# Select axes for first scatter
col1, col2 = st.columns(2)
with col1:
    x1 = st.selectbox("Select X-axis (Plot 1)", list(qc_metrics.keys()), index=0)
    y1 = st.selectbox("Select Y-axis (Plot 1)", list(qc_metrics.keys()), index=2)

with col2:
    x2 = st.selectbox("Select X-axis (Plot 2)", list(qc_metrics.keys()), index=1)
    y2 = st.selectbox("Select Y-axis (Plot 2)", list(qc_metrics.keys()), index=0)

# Create figure with 2 subplots side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1 (light red)
axes[0].scatter(
    adata.obs[qc_metrics[x1]],
    adata.obs[qc_metrics[y1]],
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
    adata.obs[qc_metrics[x2]],
    adata.obs[qc_metrics[y2]],
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
    initial_cells = adata.n_obs

    adata = adata[
        (adata.obs["n_genes_by_counts"] > min_features) &
        (adata.obs["n_genes_by_counts"] < max_features) &
        (adata.obs["pct_counts_mt"] < max_percent_mt),
        :
    ]

    # update session state so later pages see the filtered data
    st.session_state["adata"] = adata

    st.success(f"Filtered from {initial_cells} cells to {adata.n_obs} cells.")


    # Save to a temporary .h5ad file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        adata.write(tmp.name)
        tmp_path = tmp.name

    # Provide download button
    with open(tmp_path, "rb") as f:
        st.download_button(
            label="💾 Download filtered data (.h5ad)",
            data=f,
            file_name="filtered_data.h5ad",
            mime="application/octet-stream",
            key="download_filtered"  # give unique key
        )

    # Clean up temp file
    os.remove(tmp_path)