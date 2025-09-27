import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt

st.header("Step 5: Clustering and UMAP")

# --- Check if adata exists from Step 4 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 4 first.")
    st.stop()

adata = st.session_state["adata"]

# --- Subheader 1: Neighbors ---
st.subheader("Computing the neighborhood graph")

n_pcs = st.number_input(
    "Number of PCs to use for neighbors:",
    min_value=5, max_value=100, value=20, step=5
)

if st.button("Run neighbors"):
    wait_msg = st.empty()
    wait_msg.info("⏳ Please wait a moment while computing neighbors...")

    sc.pp.neighbors(adata, n_pcs=int(n_pcs))

    wait_msg.empty()
    st.success(f"✅ Nearest-neighbor graph computed (using {n_pcs} PCs).")
    st.session_state["adata"] = adata


# --- Subheader 2: Clustering + UMAP ---
st.subheader("Clustering the neighborhood graph & UMAP")

resolution = st.number_input(
    "Leiden resolution (higher = more clusters):",
    min_value=0.1, max_value=5.0, value=1.0, step=0.1, format="%.1f"
)

if st.button("Run clustering and UMAP"):
    try:
        # Placeholder for "please wait"
        wait_msg = st.empty()
        wait_msg.info("⏳ Please wait a moment while clustering and UMAP are running...")

        # Neighbors
        sc.pp.neighbors(adata, n_pcs=int(n_pcs))

        # Leiden clustering
        sc.tl.leiden(
            adata,
            resolution=float(resolution),
            random_state=0,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        adata.obs["leiden"] = adata.obs["leiden"].copy()
        adata.uns["leiden"] = adata.uns["leiden"].copy()

        # UMAP
        sc.tl.umap(adata, random_state=0)

        # Remove wait message
        wait_msg.empty()

        st.success(f"✅ Leiden clustering (resolution={resolution}, PCs={n_pcs}) and UMAP complete.")
        st.write("Number of clusters:", adata.obs["leiden"].nunique())
        st.dataframe(adata.obs["leiden"].value_counts().rename("Cell count"))

        # --- UMAP plot ---
        st.subheader("UMAP with Leiden clusters")
        sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False)
        fig = plt.gcf()
        st.pyplot(fig)
        plt.close(fig)

        # Save back to session
        st.session_state["adata"] = adata

    except ImportError as e:
        wait_msg.empty()  # clear message on error too
        st.error("Leiden requires the `igraph` and `leidenalg` packages. "
                 "Please install them via conda or pip.")
        st.exception(e)
