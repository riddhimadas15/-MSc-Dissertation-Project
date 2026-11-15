"""Rebuild tf_ra_regulations.csv (TF -> RA gene regulations) used by the RAveda app.

Input : gene_attribute_edges.txt  - ENCODE Transcription Factor Targets, downloaded from Harmonizome
        (https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets)
        columns: source (regulated gene), ..., target (transcription factor), ..., weight
Output: tf_ra_regulations.csv with columns TF, TargetGene

Only edges whose target gene is in ra_genes.csv are kept: the app filters to RA genes anyway,
and this keeps the file small. Pass --all-genes to keep every edge (~1.65 M rows).

NOTE: do not open the full edge table in Excel/LibreOffice and re-save it - spreadsheets cap at
1,048,576 rows and silently drop the rest.
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--edges", default="gene_attribute_edges.txt")
ap.add_argument("--ra-genes", default="ra_genes.csv")
ap.add_argument("--out", default="tf_ra_regulations.csv")
ap.add_argument("--all-genes", action="store_true")
a = ap.parse_args()

edges = pd.read_csv(a.edges, sep="\t", usecols=["source", "target"])
edges.columns = ["TargetGene", "TF"]
edges["TargetGene"] = edges["TargetGene"].astype(str).str.strip().str.upper()
edges["TF"] = edges["TF"].astype(str).str.strip().str.upper()

if not a.all_genes:
    ra = set(pd.read_csv(a.ra_genes)["GeneSymbol"].astype(str).str.strip().str.upper())
    edges = edges[edges["TargetGene"].isin(ra)]

edges = edges.drop_duplicates()[["TF", "TargetGene"]].sort_values(["TF", "TargetGene"])
edges.to_csv(a.out, index=False)
print(f"wrote {len(edges):,} edges, {edges.TF.nunique()} TFs -> {a.out}")
