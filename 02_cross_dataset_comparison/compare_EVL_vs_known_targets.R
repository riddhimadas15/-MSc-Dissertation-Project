# Load required libraries
library(dplyr)
library(readr)

# Read the CSV files
cve_data <- read_csv("EVL.csv")
ayur_data <- read_csv("ayur.csv")

# Use the actual column name directly (no uppercasing)
common_genes <- intersect(cve_data$gene_symbol, ayur_data$gene_symbol)

# Create a data frame from common genes
common_df <- data.frame(gene_symbol = common_genes)

# Merge with cve_data to get log2FoldChange values
common_df <- common_df %>%
  left_join(cve_data %>% select(gene_symbol, log2FoldChange), by = "gene_symbol")

# Save the result
write.csv(common_df, "common_genes_evl_ayur.csv", row.names = FALSE)

# View first few rows
head(common_df)

