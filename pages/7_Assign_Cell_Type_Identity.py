import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt

st.header("Step 7: Gene Expression and Cell Type Annotation")

# --- Check if adata exists from Step 6 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 6 first.")
    st.stop()

adata = st.session_state["adata"]

# ======================
# Part 1: Gene Expression
# ======================
st.subheader("Gene Expression on UMAP")

# Get available gene lists
all_genes = list(adata.var_names)
hvg_genes = list(adata.var[adata.var.get("highly_variable", False)].index)

gene_source = st.radio(
    "Choose gene list:",
    ["Highly variable genes", "All genes"],
    index=0
)

gene_list = hvg_genes if gene_source == "Highly variable genes" and len(hvg_genes) > 0 else all_genes

selected_genes = st.multiselect(
    "Select one or more genes:",
    options=gene_list,
    default=st.session_state.get("last_selected_genes", []),
    help="Type to search for genes of interest (e.g., CD3E, CD8A)."
)

if st.button("Plot UMAP with selected genes"):
    st.session_state["last_selected_genes"] = selected_genes
    st.session_state["show_gene_umaps"] = True

# Always replot if previously requested
if st.session_state.get("show_gene_umaps", False) and st.session_state.get("last_selected_genes"):
    colors = ["leiden"] + st.session_state["last_selected_genes"]

    # Show plots 2 per row
    for i in range(0, len(colors), 2):
        cols = st.columns(2)
        for j in range(2):
            if i + j < len(colors):
                gene = colors[i + j]
                with cols[j]:
                    sc.pl.umap(adata, color=gene, show=False, use_raw=False)
                    fig = plt.gcf()
                    st.pyplot(fig)
                    plt.close(fig)

# ======================
# Part 2: CellTypist
# ======================
st.subheader("Cell Type Annotation with CellTypist")

try:
    import celltypist
    import pandas as pd

    model_choice = st.selectbox(
        "Choose CellTypist model:",
        ["Immune_All_Low.pkl", "Immune_All_High.pkl", "Immune_Fine.pkl"]
    )

    if st.button("Run CellTypist annotation"):
        wait_placeholder = st.empty()
        wait_placeholder.info("⏳ Please wait a moment while running CellTypist...")

        # Load model
        model = celltypist.models.Model.load(model_choice)

        # Run on raw normalized data (stored in Step 3)
        prediction = celltypist.annotate(
            adata.raw.to_adata(),
            model=model,
            majority_voting=True
        )

        # --- Handle CellTypist outputs safely ---
        if hasattr(prediction, "predicted_labels"):
            labels_df = prediction.predicted_labels
            if isinstance(labels_df, pd.DataFrame):
                if "majority_voting" in labels_df.columns:
                    adata.obs["predicted_labels"] = labels_df["majority_voting"]
                else:
                    adata.obs["predicted_labels"] = labels_df.iloc[:, 0]
            else:
                adata.obs["predicted_labels"] = labels_df
        elif "predicted_labels" in prediction.adata.obs:
            adata.obs["predicted_labels"] = prediction.adata.obs["predicted_labels"]
        elif "majority_voting" in prediction.adata.obs:
            adata.obs["predicted_labels"] = prediction.adata.obs["majority_voting"]
        else:
            st.error("Could not find predicted labels in CellTypist output.")
            st.stop()

        st.success(f"✅ CellTypist annotation complete using {model_choice}")
        st.session_state["adata"] = adata

        # --- UMAP with annotation ---
        st.subheader("UMAP with predicted cell types")
        sc.pl.umap(
            adata,
            color="predicted_labels",
            legend_loc="on data",
            show=False,
            use_raw=False
        )
        fig = plt.gcf()
        st.pyplot(fig)
        plt.close(fig)

        # --- Table of predicted cell types ---
        st.subheader("Cell type counts")
        st.dataframe(adata.obs["predicted_labels"].value_counts())

        # --- Download CSV of predicted labels ---
        st.subheader("Download annotations")
        df_labels = adata.obs[["predicted_labels"]].copy()
        csv_data = df_labels.to_csv().encode("utf-8")
        st.download_button(
            label="💾 Download predicted labels (.csv)",
            data=csv_data,
            file_name="celltypist_predicted_labels.csv",
            mime="text/csv"
        )

        wait_placeholder.empty()

except ImportError:
    st.error("CellTypist is not installed. Please run: `pip install celltypist`")
