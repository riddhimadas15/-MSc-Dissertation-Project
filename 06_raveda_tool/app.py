import streamlit as st
import pandas as pd
import plotly.express as px
from pyvis.network import Network
import streamlit.components.v1 as components

# Load reference data
ayurvedic_targets = pd.read_csv("ayurvedic_targets.csv")
tf_list = pd.read_csv("transcription_factors.csv")
ra_genes = pd.read_csv("ra_genes.csv")
tf_ra_network = pd.read_csv("tf_ra_regulations.csv")

st.set_page_config(layout="wide")
st.title("🌿️ RAveda")

# Tabs
tab1, tab2 = st.tabs(["🧬 DEG Analysis", "ℹ️ About"])

# Function to display network
def show_tf_ra_network(tf_ra_df, single_tf=None):
    net = Network(height='600px', width='100%', bgcolor='#ffffff', font_color='black')

    # Add nodes and edges
    for _, row in tf_ra_df.iterrows():
        tf = row['TF']
        target = row['TargetGene']
        if single_tf and tf != single_tf:
            continue
        net.add_node(tf, label=tf, color='#FF6F61')           # TFs in red
        net.add_node(target, label=target, color='#6BAED6')   # RA genes in blue
        net.add_edge(tf, target)

    options = """
    {
      "layout": {
        "improvedLayout": true
      },
      "physics": {
        "enabled": true,
        "solver": "forceAtlas2Based",
        "forceAtlas2Based": {
          "gravitationalConstant": -50,
          "centralGravity": 0.01,
          "springLength": 100
        },
        "stabilization": {
          "enabled": true,
          "iterations": 200
        }
      }
    }
    """
    net.set_options(options)

    net.save_graph('tf_ra_network.html')

    # Stop bubble vibration after stabilization
    with open('tf_ra_network.html', 'r', encoding='utf-8') as f:
        html_content = f.read()

    html_content = html_content.replace(
        "network = new vis.Network(container, data, options);",
        """network = new vis.Network(container, data, options);
        network.once('stabilizationIterationsDone', function () {
            network.setOptions({ physics: false });
        });"""
    )

    components.html(html_content, height=650, scrolling=True)

# --- TAB 1: DEG ANALYSIS ---
with tab1:
    uploaded_file = st.file_uploader("Upload your DEG CSV file", type=['csv'])

    if uploaded_file:
        deg_df = pd.read_csv(uploaded_file)
        deg_df.columns = deg_df.columns.str.strip().str.lower()

        if 'gene' not in deg_df.columns or 'logfc' not in deg_df.columns:
            st.error("Your CSV must contain columns named 'Gene' and 'logFC'")
        else:
            # User-defined cutoff
            logfc_cutoff = st.number_input("Set logFC cutoff for upregulated genes", value=1.0)

            # Filter for upregulated genes
            upregulated_df = deg_df[deg_df['logfc'] >= logfc_cutoff]
            upregulated_genes = set(upregulated_df['gene'])

            # Match only upregulated genes to drugs
            matched_drugs = ayurvedic_targets[ayurvedic_targets['GeneTarget'].isin(upregulated_genes)]
            deg_tfs = tf_list[tf_list['GeneSymbol'].isin(upregulated_genes)]

            regulated = tf_ra_network[
                (tf_ra_network['TF'].isin(deg_tfs['GeneSymbol'])) &
                (tf_ra_network['TargetGene'].isin(ra_genes['GeneSymbol']))
            ]

            # Display Results
            st.subheader("🧪 Matched Ayurvedic Drugs (Upregulated Genes Only)")
            st.dataframe(matched_drugs.reset_index(drop=True), use_container_width=True)

            # Pie Chart: Drug -> # of GeneTargets
            st.subheader("📊 Drug Distribution by Number of Gene Targets")
            if not matched_drugs.empty:
                drug_counts = matched_drugs.groupby('DrugName')['GeneTarget'].nunique().reset_index()
                drug_counts.columns = ['DrugName', 'NumTargets']

                fig = px.pie(drug_counts, names='DrugName', values='NumTargets',
                             title="Distribution of Ayurvedic Drugs by Number of Gene Targets")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No matched drugs found for the given logFC cutoff.")

            # Unique drugs
            st.subheader("💊 Unique Ayurvedic Drugs for Treatment")
            unique_drugs = matched_drugs['DrugName'].dropna().unique()
            st.write(", ".join(sorted(unique_drugs)))

            st.subheader("🧬 Differentially Expressed TFs")
            st.dataframe(deg_tfs.reset_index(drop=True), use_container_width=True)

            st.subheader("🔗 TFs Regulating RA-Related Genes")
            st.dataframe(regulated.reset_index(drop=True), use_container_width=True)

            # Network Visualization
            if not regulated.empty:
                st.subheader("🕸️ Visualize TF–RA Gene Interaction Network")

                view_option = st.radio(
                    "Select network view:",
                    ["Show Full Network", "Show Single TF Network"]
                )

                if view_option == "Show Full Network":
                    show_tf_ra_network(regulated)

                elif view_option == "Show Single TF Network":
                    tf_options = sorted(regulated['TF'].unique())
                    selected_tf = st.selectbox("Select a Transcription Factor", tf_options)
                    tf_specific_network = regulated[regulated['TF'] == selected_tf]

                    if not tf_specific_network.empty:
                        show_tf_ra_network(tf_specific_network, single_tf=selected_tf)
                    else:
                        st.info("No RA-related target genes found for the selected TF.")
            else:
                st.info("No TF–RA interactions found in your DEGs.")

# --- TAB 2: ABOUT ---
with tab2:
    st.markdown("### 🔍 About This Tool")
    st.markdown("""
    This interactive tool helps researchers explore potential Ayurvedic treatments for Rheumatoid Arthritis (RA) 
    based on gene expression data. It performs the following:

    - Matches differentially expressed genes (DEGs) with gene targets of Ayurvedic compounds.
    - Filters for upregulated genes based on a logFC threshold.
    - Identifies transcription factors (TFs) among DEGs.
    - Checks which of those TFs regulate RA-related genes.
    - Visualizes TF–RA gene interactions either fully or for selected TFs using a network graph.
    - Shows a pie chart of Ayurvedic drugs and how many upregulated genes they target.

    **Input:** CSV file with columns `Gene` and `logFC` (case-insensitive).

    **Developed by Riddhima – April 2025**
    """)

