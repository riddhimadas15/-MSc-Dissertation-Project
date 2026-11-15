# 3. Transcription-factor network analysis

Asks which DEGs are themselves transcription factors (TFs), which genes they regulate, and whether those targets are known RA genes.

**Data sources**
- TF → target interactions: ENCODE Transcription Factor Targets via [Harmonizome](https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets) (`gene_attribute_edges.txt`, ~65 MB, **download separately** into this folder; `attribute_list_entries.txt` = the TF list is included; see `explanation_of_harmonizome_files.txt`).
- `all_RA_genes_source.csv`: 212 RA-associated genes curated from KEGG and the literature, with the source of each.
- DEGs: `GSM265020_filtered_genes_1.csv` (GSE10500, from stage 1).

**Run order** (from this folder)

| Script | Purpose | Main output |
|---|---|---|
| `01_find_DEG_TFs_and_targets.R` | DEGs that are TFs, their targets, and TFs regulating the non-TF DEGs | `results/DEG_TF_Targets.csv`, `results/DEG_Regulating_TFs.csv` |
| `02_overlap_TF_targets_with_RA_genes.R` | Overlap of each TF's targets with the RA gene list | `common_genes*.csv` (thesis Table 4.3) |
| `03_plot_TF_gene_network.R` | Directed TF–gene network (igraph) | plot |
| `04_plot_single_TF_network.R` | Single-TF regulatory network, e.g. SPI1 (ggraph) | plot |
| `05_heatmap_common_RA_genes.R` | Heatmap of expression of DEGs that are also RA genes | `expression_plot_66.png` |
| `06_overlap_DEGs_with_RA_genes.R` | DEGs ∩ RA gene list at different fold-change cut-offs | `common_genes_*.csv` |

`results/DEG_TF_Gene_Connections.csv` is the per-TF summary used by script 02. Shipped result tables and a heatmap are in `results/`.
