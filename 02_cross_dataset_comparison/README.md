# 2. Cross-dataset comparison

Checks which genes are dysregulated consistently across datasets, and which of them are known Ayurvedic-bioactive targets.

- Microarray DEG lists (GSE10500, GSE93272, GSE55235, GSE206848, GSE77298) were compared gene by gene; genes found in at least three datasets were kept for network analysis (thesis §3.3).
- RNA-seq (GSE89408) contrasts: **CVE** = control vs early RA, **CVL** = control vs late RA, **EVL** = early vs late RA.

## Scripts
Each script reads two tables with a `gene_symbol` column, intersects them with `dplyr`, and writes the shared genes (optionally with log2 fold change):

| Script | Compares |
|---|---|
| `compare_CVE_vs_known_targets.R` | RNA-seq CVE DEGs ↔ known bioactive-target genes |
| `compare_CVL_vs_known_targets.R` | RNA-seq CVL DEGs ↔ known bioactive-target genes |
| `compare_EVL_vs_known_targets.R` | RNA-seq EVL DEGs ↔ known bioactive-target genes |
| `compare_CVE_vs_CVL.R` | CVE ↔ CVL DEGs |

## Inputs not included
The RNA-seq DESeq2 result tables (`CVE.csv`, `CVL.csv`, `EVL.csv`) and the other datasets' DEG tables are not redistributed here. To reuse the scripts, supply your own tables and edit the file names at the top of each script. The known-target list is in [`../04_bioactive_target_mapping/data/ayurvedic_targets.csv`](../04_bioactive_target_mapping/data/ayurvedic_targets.csv).
