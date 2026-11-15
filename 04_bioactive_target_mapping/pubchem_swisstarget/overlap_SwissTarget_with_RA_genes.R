# Load required libraries
library(dplyr)
library(readr)

# Read the CSV files
cve_data <- read_csv("genes_RA.csv")
ayur_data <- read_csv("SwissTargetPrediction.csv")

# Use the actual column name directly (no uppercasing)
common_genes <- intersect(cve_data$gene_symbol, ayur_data$Common_name)

# Create a data frame from common genes
common_df <- data.frame(gene_symbol = common_genes)

# Save the result
write.csv(common_df, "common_genes_cve_ayur.csv", row.names = FALSE)

# View first few rows
head(common_df)




