import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

st.header("Step 6: Differential Expression (Marker Genes)")

# --- Check if adata exists from Step 5 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 5 first.")
    st.stop()

adata = st.session_state["adata"]

if "leiden" not in adata.obs:
    st.error("No clustering found. Please run Step 5 first.")
    st.stop()

# --- Select comparison mode ---
mode = st.radio(
    "Choose comparison mode:",
    ["Cluster vs all other clusters", "Cluster vs cluster"]
)

clusters = sorted(adata.obs["leiden"].unique())

run_deg = False
cluster = None
cluster1 = None
cluster2 = None

if mode == "Cluster vs all other clusters":
    cluster = st.selectbox("Select cluster:", clusters)
    if st.button("Run DE analysis"):
        run_deg = True

elif mode == "Cluster vs cluster":
    cluster1 = st.selectbox("Select cluster 1:", clusters, index=0)
    cluster2 = st.selectbox("Select cluster 2:", clusters, index=1)
    if st.button("Run DE analysis"):
        run_deg = True

# --- Run DE ---
if run_deg:
    wait_msg = st.empty()
    wait_msg.info("⏳ Please wait a moment while differential expression analysis is running...")

    if mode == "Cluster vs all other clusters":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster],
            reference="rest",
            method="wilcoxon"
        )
        st.success(f"✅ Marker genes for cluster {cluster} vs all others computed.")

    elif mode == "Cluster vs cluster":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster1],
            reference=cluster2,
            method="wilcoxon"
        )
        st.success(f"✅ Marker genes for cluster {cluster1} vs cluster {cluster2} computed.")

    wait_msg.empty()  # remove "please wait" message
    st.session_state["adata"] = adata

# --- Show results if available ---
if "rank_genes_groups" in adata.uns:
    st.subheader("Top marker genes (Scanpy plot)")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # Convert results to DataFrame
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "names": result["names"][g],
            "scores": result["scores"][g],
            "logfoldchanges": result["logfoldchanges"][g],
            "pvals": result["pvals"][g],
            "pvals_adj": result["pvals_adj"][g],
        })
        df["cluster"] = g
        dfs.append(df)
    df_out = pd.concat(dfs)

    # Download button
    st.download_button(
        label="💾 Download DE results (.csv)",
        data=df_out.to_csv(index=False).encode("utf-8"),
        file_name="DE_results.csv",
        mime="text/csv"
    )
