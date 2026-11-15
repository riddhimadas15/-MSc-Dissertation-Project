# Install required packages if not already installed
if (!require(igraph)) install.packages("igraph", dependencies=TRUE)
if (!require(ggraph)) install.packages("ggraph", dependencies=TRUE)
if (!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)

# Load libraries
library(igraph)
library(ggraph)
library(ggplot2)

# Load CSV file containing SPI1 and its target genes
file_path <- "results/DEG_TF_Targets.csv"  # Ensure this file exists in the working directory
edges <- read.csv(file_path, stringsAsFactors = FALSE)

# Ensure column names match your dataset
colnames(edges) <- c("TF", "Target_Gene")  # Update column names if needed

# Check if data is loaded correctly
if (nrow(edges) == 0) {
  stop("Error: CSV file is empty or incorrectly formatted.")
}

# Convert TF and Target_Gene to character type
edges$TF <- as.character(edges$TF)
edges$Target_Gene <- as.character(edges$Target_Gene)

# Create graph from CSV data
gene_network <- graph_from_data_frame(edges, directed = TRUE)

# Ensure graph is not empty
if (vcount(gene_network) == 0 || ecount(gene_network) == 0) {
  stop("Error: Graph is empty. Check input data.")
}

# Define RA-associated genes (Example list, update as needed)
RA_genes <- c("VEGFA", "TNFAIP3", "S100A12", "CD14", "F13A1")  

# Assign colors based on node type
V(gene_network)$color <- ifelse(V(gene_network)$name == "SPI1", "purple",  # TF SPI1 in purple
                         ifelse(V(gene_network)$name %in% RA_genes, "red", "blue"))  # RA genes in red, others in blue

# Assign node sizes based on degree (higher degree = bigger node)
V(gene_network)$size <- degree(gene_network) * 1.5  

# Use a circular layout for a cleaner structure
p <- ggraph(gene_network, layout = "circle") +  
  geom_edge_link(aes(alpha = 0.5), edge_colour = "gray", width = 0.4) +  # Faded, thinner edges
  geom_node_point(aes(size = size, color = color)) +  # Node size & color based on degree
  geom_node_text(aes(label = name), repel = TRUE, size = 4, vjust = 1.5) +  # Avoid overlapping text
  scale_color_identity() +  # Use predefined colors
  theme_void() +
  ggtitle("SPI1 Regulatory Network with RA-Associated Genes")

# Print plot
print(p)

