import csv
import re

def parse_herb_drug_data(txt_file, output_csv):
    with open(txt_file, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]

    data = []
    current_herb = None
    current_drug = None
    gene_pattern = re.compile(r"^[A-Z0-9]+$")

    for line in lines:
        # Herb name (bold in original, assume lines with no tabs or parenthesis and multiple words)
        if re.match(r"^[A-Za-z\s/]+$", line) and not '(' in line:
            current_herb = line.strip()
            continue

        # Drug + SMILES line (tab-separated or multiple spaces)
        if re.search(r'\s{2,}|\t', line):
            parts = re.split(r'\s{2,}|\t+', line)
            if len(parts) >= 2:
                drug_name_raw = parts[0].strip()
                current_drug = drug_name_raw.split('(')[0].strip()
                continue

        # Gene target line
        if current_drug and gene_pattern.match(line):
            data.append([current_drug, line])

    # Write to CSV
    with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['DrugName', 'GeneTarget'])
        writer.writerows(data)

    print(f"✅ Extraction complete! File saved as: {output_csv} with {len(data)} rows.")

# Example usage:
parse_herb_drug_data('data/herb_bioactive_target_notes.txt', 'data/drug_gene_targets.csv')

