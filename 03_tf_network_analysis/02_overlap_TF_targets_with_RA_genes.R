library(dplyr)

# Load the CSV files
filtered_genes <- read.csv("results/DEG_TF_Gene_Connections.csv", stringsAsFactors = FALSE)
all_ra_genes <- read.csv("all_RA_genes_source.csv", stringsAsFactors = FALSE)

# Ensure the column names match (adjust if needed)
gene_column_name <- "Genes"  # Change this if your column name is different

# Function to detect common genes
detect_common_genes <- function(genes_column, all_ra_genes) {
  common_genes <- intersect(genes_column, all_ra_genes[["RA_genes"]])
  return(common_genes)
}

# Apply function to check for each TF's gene list and add Num_Connections & Common_new columns
filtered_genes <- filtered_genes %>% 
  mutate(
    Common_Genes = sapply(strsplit(Genes, ", "), function(genes) {
      paste(detect_common_genes(genes, all_ra_genes), collapse = ", ")
    }),
    Num_Connections = sapply(strsplit(Genes, ", "), length),
    Common_new = sapply(strsplit(Genes, ", "), function(genes) {
      length(detect_common_genes(genes, all_ra_genes))
    })
  )

# Save to a new CSV file
write.csv(filtered_genes, "common_genes.csv", row.names = FALSE)

cat("Analysis completed. Common genes saved in common_genes.csv.\n")


# Load the CSV files
filtered_genes <- read.csv("results/DEG_TF_Gene_Connections_70.csv", stringsAsFactors = FALSE)
all_ra_genes <- read.csv("all_RA_genes_source.csv", stringsAsFactors = FALSE)

# Ensure the column names match (adjust if needed)
gene_column_name <- "Genes"  # Change this if your column name is different

# Function to detect common genes
detect_common_genes <- function(genes_column, all_ra_genes) {
  common_genes <- intersect(genes_column, all_ra_genes[["RA_genes"]])
  return(common_genes)
}

# Apply function to check for each TF's gene list and add Num_Connections & Common_new columns
filtered_genes <- filtered_genes %>% 
  mutate(
    Common_Genes = sapply(strsplit(Genes, ", "), function(genes) {
      paste(detect_common_genes(genes, all_ra_genes), collapse = ", ")
    }),
    Num_Connections = sapply(strsplit(Genes, ", "), length),
    Common_new = sapply(strsplit(Genes, ", "), function(genes) {
      length(detect_common_genes(genes, all_ra_genes))
    })
  )

# Save to a new CSV file
write.csv(filtered_genes, "common_genes_1.5.csv", row.names = FALSE)

cat("Analysis completed. Common genes saved in common_genes.csv.\n")


