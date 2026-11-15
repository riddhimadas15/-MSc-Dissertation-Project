library(dplyr)

# Load the CSV files for logfc 1
filtered_genes <- read.csv("GSM265020_filtered_genes_1.csv", stringsAsFactors = FALSE)
all_ra_genes <- read.csv("all_RA_genes_source.csv", stringsAsFactors = FALSE)

# Ensure the column name matches (adjust if needed)
gene_column_name <- "SYMBOL"  # Change this if your column name is different
gene_column_name_1 <- "RA_genes"

# Find common genes
common_genes_1 <- intersect(filtered_genes[[gene_column_name]], all_ra_genes[[gene_column_name_1]])

# Save to a new CSV file
common_genes_df <- data.frame(Gene = common_genes_1)
write.csv(common_genes_df, "common_genes_1.csv", row.names = FALSE)

cat("Found", length(common_genes_1), "common genes. Saved to common_genes.csv.\n")

# Load the CSV files _ fpr logfc1.5
filtered_genes1 <- read.csv("GSM265020_filtered_genes_1.5.csv", stringsAsFactors = FALSE)
all_ra_genes1 <- read.csv("all_RA_genes_source.csv", stringsAsFactors = FALSE)

# Ensure the column name matches (adjust if needed)
gene_column_name2 <- "SYMBOL"  # Change this if your column name is different
gene_column_name_2 <- "RA_genes"

# Find common genes
common_genes1 <- intersect(filtered_genes1[[gene_column_name2]], all_ra_genes1[[gene_column_name_2]])

# Save to a new CSV file
common_genes_df_1 <- data.frame(Gene = common_genes1)
write.csv(common_genes_df_1, "common_genes1.5.csv", row.names = FALSE)

cat("Found", length(common_genes1), "common genes. Saved to common_genes1.5.csv.\n")


# Load the CSV files _ fpr logfc1.7
filtered_genes_1.7 <- read.csv("GSM265020_filtered_genes_1.7.csv", stringsAsFactors = FALSE)
all_ra_genes1.7 <- read.csv("all_RA_genes_source.csv", stringsAsFactors = FALSE)

# Ensure the column name matches (adjust if needed)
gene_column_name1.7 <- "SYMBOL"  # Change this if your column name is different
gene_column_name_1.7<- "RA_genes"

# Find common genes
common_genes2 <- intersect(filtered_genes_1.7[[gene_column_name1.7]], all_ra_genes1.7[[gene_column_name_1.7]])

# Save to a new CSV file
common_genes_df_2 <- data.frame(Gene = common_genes2)
write.csv(common_genes_df_2, "common_genes1.7.csv", row.names = FALSE)

cat("Found", length(common_genes2), "common genes. Saved to common_genes1.7.csv.\n")


# Load the CSV files _ fpr logfc0.76 aka linear 1.7
filtered_genes_linear_1.7 <- read.csv("GSM265020_filtered_genes_linear_1.7", stringsAsFactors = FALSE)
all_ra_genes_linear1.7 <- read.csv("all_RA_genes_source.csv", stringsAsFactors = FALSE)

# Ensure the column name matches (adjust if needed)
gene_column_name_linear_1.7 <- "SYMBOL"  # Change this if your column name is different
gene_column_name_line_1.7<- "RA_genes"

# Find common genes
common_genes_lin_1.7 <- intersect(filtered_genes_linear_1.7[[gene_column_name_linear_1.7]], all_ra_genes_linear1.7[[gene_column_name_line_1.7]])

# Save to a new CSV file
common_genes_df_lin_1.7 <- data.frame(Gene = common_genes_lin_1.7)
write.csv(common_genes_df_lin_1.7, "common_genes_linear_1.7.csv", row.names = FALSE)

cat("Found", length(common_genes_lin_1.7), "common genes. Saved to common_genes_lin_1.7.csv.\n")




