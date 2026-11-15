# Load required package
library(igraph)

# Load DEGs from limma results
deg_data <- read.csv("GSM265020_filtered_genes_1.csv", header = TRUE, stringsAsFactors = FALSE)

# Load transcription factor data
tf_data <- read.table("attribute_list_entries.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Convert GeneSym column to uppercase for case-insensitive matching
tf_data$GeneSym <- toupper(tf_data$GeneSym)

# Load gene-TF interaction network data
gene_tf_network_data <- read.table("gene_attribute_edges.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Ensure correct column names
colnames(gene_tf_network_data) <- c("source", "source_desc", "source_id", "target", "target_desc", "target_id", "weight")

# Keep only relevant columns
gene_tf_network_data <- gene_tf_network_data[, c("source", "target")]

# Convert all gene names to uppercase for consistency
gene_tf_network_data$source <- toupper(gene_tf_network_data$source)
gene_tf_network_data$target <- toupper(gene_tf_network_data$target)
deg_data$SYMBOL <- toupper(deg_data$SYMBOL)

# Remove missing values from DEGs
deg_genes <- na.omit(deg_data$SYMBOL)

# Create graph from gene-TF interaction table
gene_tf_graph <- graph_from_data_frame(d = gene_tf_network_data, directed = TRUE)

# Ensure graph vertices have valid names
V(gene_tf_graph)$name <- as.character(V(gene_tf_graph)$name)

# Identify which DEGs are transcription factors
deg_tf_list <- intersect(deg_genes, tf_data$GeneSym)  # Compare with TF list

# Print debug information
print("Identified DEG Transcription Factors:")
print(deg_tf_list)

# Identify which genes they regulate
tf_to_gene_list <- data.frame(TF = character(), Target_Gene = character(), stringsAsFactors = FALSE)
for (tf in deg_tf_list) {
  if (tf %in% V(gene_tf_graph)$name) {  
    targets <- names(neighbors(gene_tf_graph, tf, mode = "in"))  # Compare with target column
    if (length(targets) > 0) {
      tf_to_gene_list <- rbind(tf_to_gene_list, data.frame(TF = rep(tf, length(targets)), Target_Gene = targets))
    }
  } else {
    print(paste("TF", tf, "not found in the network graph"))
  }
}

# Identify non-TF DEGs and their regulators
non_tf_degs <- setdiff(deg_genes, tf_data$GeneSym)  # Compare with TF list

gene_to_tf_list <- data.frame(Gene = character(), Regulating_TF = character(), stringsAsFactors = FALSE)
for (gene in non_tf_degs) {
  if (gene %in% V(gene_tf_graph)$name) {  
    regulators <- names(neighbors(gene_tf_graph, gene, mode = "out"))  # Compare with source column
    if (length(regulators) > 0) {
      gene_to_tf_list <- rbind(gene_to_tf_list, data.frame(Gene = rep(gene, length(regulators)), Regulating_TF = regulators))
    }
  } else {
    print(paste("Gene", gene, "not found in the network graph"))
  }
}

# Print results
print("Transcription Factors in DEGs and their target genes:")
print(tf_to_gene_list)

print("Non-TF DEGs and the TFs regulating them:")
print(gene_to_tf_list)

# Save results only if they are non-empty
if (nrow(tf_to_gene_list) > 0) {
  write.csv(tf_to_gene_list, "results/DEG_TF_Targets.csv", row.names = FALSE)
}

if (nrow(gene_to_tf_list) > 0) {
  write.csv(gene_to_tf_list, "results/DEG_Regulating_TFs.csv", row.names = FALSE)
}
