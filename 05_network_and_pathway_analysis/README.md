# 5. Network and pathway analysis

- **`centrality_analysis.py`**: builds a graph from a STRING network export and ranks genes by degree, betweenness and closeness centrality; saves a table and histograms. Used for (a) the bioactives' target network and (b) the network of consensus DEGs.
  ```bash
  python centrality_analysis.py --input data/string_interactions_bioactive_targets.tsv --name bioactive_targets
  ```
  In the included network, *TNF* has the highest degree, followed by *IL10*, *CXCL8*, *IL1B*, *CXCL1* and *CCL5*.
- **`common_pathways.py`**: intersects the KEGG pathways enriched for the common DEGs (`data/enrichment_KEGG.tsv`) with those enriched for the bioactives' targets (`data/enrichment_KEGG_targets.tsv`). 64 pathways are shared, including Th17 differentiation, TNF, NF-κB, MAPK and IL-17 signalling.
- **`formulation_networks/`**: formulation → bioactive → target and formulation → pathway graphs for Ashwashila, Articulin F and Triphala Guggul (`networkx` + `matplotlib`).
- `results/figures/`: betweenness and closeness centrality plots used in the thesis.

STRING settings: highest confidence (0.9), first-shell interactors only. Enrichment was run in DAVID (KEGG, BH-adjusted p < 0.05).
