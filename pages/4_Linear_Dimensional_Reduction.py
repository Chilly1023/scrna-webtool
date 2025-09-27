import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import math

st.header("Step 4: Linear Dimensional Reduction (PCA)")

# --- Check if adata exists from Step 3 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 3 first.")
    st.stop()

adata = st.session_state["adata"]

# --- Run PCA ---
if st.button("Run PCA"):
    wait_placeholder = st.empty()
    wait_placeholder.info("⏳ Please wait a moment while PCA is running...")

    sc.tl.pca(adata, svd_solver="arpack")

    st.success("✅ PCA computed and stored in AnnData.")
    st.session_state["adata"] = adata

    # remove the "please wait" message
    wait_placeholder.empty()

# --- If PCA done, show plots ---
if "X_pca" in adata.obsm_keys():
    # ---- Elbow plot ----
    st.subheader("Variance explained (Elbow plot)")
    n_pcs_elbow = st.number_input(
        "Number of PCs to display in elbow plot:",
        min_value=5, max_value=100, value=20, step=5
    )
    sc.pl.pca_variance_ratio(adata, log=False, n_pcs=n_pcs_elbow, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # ---- PCA Loadings in pairs ----
    st.subheader("Top contributing genes per PC (pairs)")
    n_pcs_loadings = st.number_input(
        "Number of PCs to plot loadings for:",
        min_value=2, max_value=20, value=6, step=2  # must be even
    )

    n_pairs = math.ceil(n_pcs_loadings / 2)
    for i in range(n_pairs):
        pc1 = 2 * i + 1
        pc2 = 2 * i + 2
        if pc2 <= n_pcs_loadings:
            sc.pl.pca_loadings(
                adata,
                components=(pc1, pc2),
                include_lowest=True,
                show=False
            )
            fig = plt.gcf()
            st.pyplot(fig)
            plt.close(fig)

else:
    st.info("👉 Run PCA first to view elbow plot and PC loadings.")
