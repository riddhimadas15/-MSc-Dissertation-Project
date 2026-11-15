import pandas as pd

# Replace with your actual TSV file paths
file_paths = ['data/enrichment_KEGG.tsv', 'data/enrichment_KEGG_targets.tsv']

# Read pathway columns into sets
pathway_sets = []
for path in file_paths:
    df = pd.read_csv(path, sep='\t')  # TSV files use tab as separator
    pathways = set(df['pathway'].dropna().astype(str).str.strip())  # Convert to string, remove whitespace
    pathway_sets.append(pathways)

# Find intersection of all sets
common_pathways = set.intersection(*pathway_sets)

# Print or save the result
print(f"Number of common pathways: {len(common_pathways)}")
print("Common pathways:", common_pathways)

# Optional: Save to a new TSV
pd.DataFrame(sorted(common_pathways), columns=["common_pathways"]).to_csv("results/common_pathways.tsv", sep='\t', index=False)

