import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

st.title("🧭 Gene Expression and Cell Type Annotation")

st.markdown("""
In this step, we explore **gene expression patterns** on UMAP and annotate clusters with **cell type labels**.

The workflow includes:
1. **Gene expression on UMAP** – visualize expression of marker or user-selected genes.  
2. **Automatic marker gene detection** – compute top genes distinguishing clusters.  
3. **Cell type annotation (CellTypist)** – automatically predict cell types using pretrained models.  

👉 This step helps us assign **biological meaning** to the clusters identified earlier.
""")

# --- Check if adata exists from Step 6 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 6 first.")
    st.stop()

adata = st.session_state["adata"]

# =========================================================
# Part 1: Gene Expression Visualization
# =========================================================
st.subheader("📌 Step 1: Gene Expression Visualization")

gene_source = st.radio(
    "Choose gene list:",
    ["Highly variable genes", "All genes"],
    index=0
)

# Get available gene lists
all_genes = adata.var_names.tolist()


if "highly_variable" in adata.var.columns:
    hvg_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
elif "highly_variable_nbatches" in adata.var.columns:  
    hvg_mask = (adata.var["highly_variable_nbatches"] > 0).to_numpy(dtype=bool)
else:
    hvg_mask = np.zeros(adata.n_vars, dtype=bool)

hvg_genes = adata.var_names[hvg_mask].tolist()


gene_list = hvg_genes if (gene_source == "Highly variable genes" and len(hvg_genes) > 0) else all_genes
if gene_source == "Highly variable genes" and len(hvg_genes) == 0:
    st.warning("No HVG mask found on this AnnData; falling back to all genes.")




gene_list = hvg_genes if gene_source == "Highly variable genes" and len(hvg_genes) > 0 else all_genes


st.info("""
💡 **Tip:** Marker genes highlight specific cell types. Common examples in PBMC:  
- **CST3** → dendritic cell / monocyte marker  
- **NKG7** → NK / cytotoxic T cell marker  
- **MS4A1** → B cell marker  
- **CD3D** → T cell marker  
- **PPBP** → Platelet marker  
""")


marker_genes = ["CST3", "NKG7", "PPBP"]


selected_genes = st.multiselect(
    "Select one or more genes:",
    options=adata.var_names.tolist(),
    default=marker_genes,
    help="Choose genes to visualize"
)

st.markdown("""
You can explore how genes are expressed across clusters using different visualization methods:  
- **UMAP** → shows spatial patterns of gene expression in 2D.  
- **Violin plots** → shows the **distribution** of expression across clusters.  
""")

plot_type = st.radio(
    "Choose visualization type:",
    ["UMAP", "Violin plots" ] #, "Marker gene table"]
)

if st.button("Generate plots"):
    if plot_type == "UMAP":
        for gene in selected_genes:
            st.subheader(f"UMAP: {gene}")
            sc.pl.umap(adata, color=gene, show=False, use_raw=False)
            fig = plt.gcf()
            st.pyplot(fig)
            plt.close(fig)

    elif plot_type == "Violin plots":
        for gene in selected_genes:
            st.subheader(f"Violin plot: {gene}")
            sc.pl.violin(adata, keys=gene, groupby="leiden", show=False)
            fig = plt.gcf()
            st.pyplot(fig)
            plt.close(fig)

    # elif plot_type == "Marker gene table":
    #     sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    #     result = adata.uns["rank_genes_groups"]
    #     groups = result["names"].dtype.names
    #     dfs = []
    #     for g in groups:
    #         df = pd.DataFrame({
    #             "names": result["names"][g][:10],
    #             "scores": result["scores"][g][:10],
    #             "logfoldchanges": result["logfoldchanges"][g][:10],
    #             "pvals_adj": [f"{x:.2e}" for x in result["pvals_adj"][g][:10]]
    #         })
    #         df["cluster"] = g
    #         dfs.append(df)
    #     df_out = pd.concat(dfs)
    #     st.dataframe(df_out)

# =========================================================
# Part 2: Automatic Marker Gene Detection
# =========================================================
st.subheader("📌 Step 2: Automatic Marker Gene Detection")

st.markdown("""
Here we use **differential expression analysis** to automatically find **marker genes** for each cluster.  
These are genes that are **highly expressed in one cluster compared to others**.
""")

top_n = st.slider("Number of top marker genes per cluster:", min_value=3, max_value=20, value=5)

if st.button("Find marker genes"):
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )
    st.success(f"✅ Computed marker genes for all clusters (top {top_n}).")
    st.session_state["adata"] = adata

    # Plot Scanpy result
    sc.pl.rank_genes_groups(adata, n_genes=top_n, sharey=False, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # Convert to DataFrame
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "names": result["names"][g][:top_n],
            "scores": result["scores"][g][:top_n],
            "logfoldchanges": result["logfoldchanges"][g][:top_n],
            "pvals_adj": [f"{x:.2e}" for x in result["pvals_adj"][g][:top_n]]
        })
        df["cluster"] = g
        dfs.append(df)
    df_out = pd.concat(dfs)

    st.dataframe(df_out)

    st.download_button(
        label="💾 Download marker genes (.csv)",
        data=df_out.to_csv(index=False).encode("utf-8"),
        file_name="marker_genes.csv",
        mime="text/csv"
    )

# =========================================================
# Part 3: Cell Type Annotation (CellTypist)
# =========================================================
st.subheader("📌 Step 3: Cell Type Annotation (CellTypist)")

st.markdown("""
So far, clustering has only given us **cluster numbers** (0, 1, 2 …).  
These numbers show groups of similar cells, but they don’t tell us **what type of cells** they are.

Here we use **CellTypist**, a machine learning tool trained on thousands of annotated single-cell datasets,  
to automatically predict the **biological identity** of each cell (e.g., T cells, B cells, NK cells, monocytes, platelets).

👉 In short:  
- **Leiden clustering** = mathematical grouping of similar cells.  
- **CellTypist** = translate those groups into known **cell types**.  

The results will be shown on the UMAP plot and as a summary table.
""")
