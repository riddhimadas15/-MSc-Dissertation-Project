# 6. RAveda: DEG → Ayurvedic drug / TF lookup

A small Streamlit app that packages parts of stages 3 and 4 for reuse on any RA differential-expression result.

**Input:** a CSV with columns `Gene` and `logFC` (see `example_input_DEGs.csv`).
**Output:** upregulated genes (user-set logFC cut-off) matched to Ayurvedic bioactives and their target counts (pie chart); DEG transcription factors; TFs regulating RA genes; interactive TF → RA-gene network (full or per TF).

```bash
pip install -r ../requirements.txt
streamlit run app.py
```

| File | Content |
|---|---|
| `ayurvedic_targets.csv` | bioactive → gene-target table (stage 4) |
| `transcription_factors.csv` | list of TFs (ENCODE) |
| `ra_genes.csv` | curated RA gene list |
| `tf_ra_regulations.csv` | TF → RA-gene regulations (rebuilt with `build_tf_ra_regulations.py`) |

To regenerate the regulation table, download `gene_attribute_edges.txt` from Harmonizome (see stage 3) and run `python build_tf_ra_regulations.py --edges /path/to/gene_attribute_edges.txt`.
