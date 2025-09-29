import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

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
# Part 1: Gene Expression on UMAP
# =========================================================
st.subheader("📌 Step 1: Gene Expression on UMAP")

# Marker gene tip
st.info("""
💡 **Tip:** Marker genes are genes whose expression highlights specific cell types.  
Here are some commonly used marker genes in PBMC data:

- **CST3** → dendritic cell / monocyte marker  
- **NKG7** → NK cell / cytotoxic T cell marker  
- **MS4A1** → B cell marker  
- **CD3D** → T cell marker  
- **PPBP** → Platelet marker  
""")

# Default marker list
marker_genes = ["CST3", "NKG7", "PPBP"]

# Gene list source
gene_source = st.radio(
    "Choose gene list:",
    ["Highly variable genes", "All genes"],
    index=0,
    help="""
- **Highly variable genes (HVGs)**: Focus only on the most informative genes (faster, more concise list).  
- **All genes**: Full list of genes in the dataset. Choose this if your marker gene is not in HVGs.  
"""
)

# Build gene list
all_genes = list(adata.var_names)
hvg_genes = list(adata.var[adata.var.get("highly_variable", False)].index)
gene_list = hvg_genes if gene_source == "Highly variable genes" and len(hvg_genes) > 0 else all_genes

# Select genes
selected_genes = st.multiselect(
    "Select one or more genes to visualize:",
    options=gene_list,
    default=marker_genes,
    help="Choose from the list of genes to color cells on UMAP."
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


# if st.button("Plot UMAP with selected genes"):
#     st.session_state["last_selected_genes"] = selected_genes
#     st.session_state["show_gene_umaps"] = True

# # Re-plot if previously requested
# if st.session_state.get("show_gene_umaps", False) and st.session_state.get("last_selected_genes"):
#     for gene in st.session_state["last_selected_genes"]:
#         sc.pl.umap(adata, color=gene, show=False, use_raw=False)
#         fig = plt.gcf()
#         st.pyplot(fig)
#         plt.close(fig)

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
