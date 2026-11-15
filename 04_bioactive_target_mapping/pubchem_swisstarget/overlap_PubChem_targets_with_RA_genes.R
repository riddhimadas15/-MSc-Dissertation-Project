# Load required libraries
library(dplyr)
library(readr)

# Read the CSV files
cve_data <- read_csv("genes_RA.csv")
ayur_data <- read_csv("p1.csv")

# Filter ayur_data to keep only Homo sapiens rows
ayur_filtered <- ayur_data %>% filter(Taxonomy == "Homo sapiens (human)")

# Find common gene symbols
common_genes <- intersect(cve_data$gene_symbol, ayur_filtered$Gene)

# Create a data frame from common genes
common_df <- data.frame(gene_symbol = common_genes)

# Save the result
write.csv(common_df, "common_genes_cve_pub.csv", row.names = FALSE)

# View first few rows
head(common_df)

