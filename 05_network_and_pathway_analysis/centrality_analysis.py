"""Hub-gene analysis of a STRING protein-protein interaction network.

Computes degree, betweenness and closeness centrality for every node, writes a ranked table,
and saves histograms of the betweenness / closeness distributions.

Used twice in the thesis (Sections 4.2.2 and 4.2.3):
  * network of the bioactives' target genes     -> data/string_interactions_bioactive_targets.tsv
  * network of the common DEGs across datasets  -> export the STRING 'network edges' TSV and pass it with --input

Usage:
    python centrality_analysis.py --input data/string_interactions_bioactive_targets.tsv --name bioactive_targets
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--input", required=True, help="STRING TSV export with Protein_1 / Protein_2 columns")
ap.add_argument("--name", default="network", help="prefix for output files")
ap.add_argument("--outdir", default="results")
ap.add_argument("--top", type=int, default=10)
a = ap.parse_args()

out = Path(a.outdir)
(out / "figures").mkdir(parents=True, exist_ok=True)

df = pd.read_csv(a.input, sep="\t")
G = nx.from_pandas_edgelist(df, source="Protein_1", target="Protein_2")

table = pd.DataFrame({
    "degree_centrality": nx.degree_centrality(G),
    "betweenness_centrality": nx.betweenness_centrality(G),
    "closeness_centrality": nx.closeness_centrality(G),
})
table["degree"] = pd.Series(dict(G.degree()))
table = table.sort_values("degree_centrality", ascending=False)
table.index.name = "gene"
table.to_csv(out / f"{a.name}_centralities.csv")

print(f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
print(f"Top {a.top} hubs by degree centrality:")
print(table.head(a.top)[["degree", "degree_centrality"]].round(4).to_string())

for col, color in [("betweenness_centrality", "skyblue"), ("closeness_centrality", "lightgreen")]:
    plt.figure(figsize=(10, 5))
    plt.hist(table[col], bins=30, color=color, edgecolor="black")
    plt.title(f"Histogram of {col.replace('_', ' ').title()} ({a.name})")
    plt.xlabel(col.replace("_", " ").title()); plt.ylabel("Frequency"); plt.grid(True); plt.tight_layout()
    plt.savefig(out / "figures" / f"{a.name}_{col}_hist.png", dpi=200); plt.close()
