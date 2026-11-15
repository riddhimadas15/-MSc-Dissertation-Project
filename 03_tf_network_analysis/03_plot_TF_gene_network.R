# Load necessary library
library(igraph)

# Read data files
df_genes <- read.csv("GSM265020_filtered_genes_1.csv", header = TRUE)
tf_data <- read.table("attribute_list_entries.txt", sep = "\t", header = TRUE)
gene_tf_network_data <- read.table("gene_attribute_edges.txt", sep = "\t", header = TRUE)

# Convert all names to uppercase to prevent mismatch errors
df_genes$SYMBOL <- toupper(df_genes$SYMBOL)
tf_data$GeneSym <- toupper(tf_data$GeneSym)
gene_tf_network_data$source <- toupper(gene_tf_network_data$source)
gene_tf_network_data$target <- toupper(gene_tf_network_data$target)

# Create a vector of genes (from the gene list)
genes_only <- df_genes$SYMBOL

# Create a vector of transcription factors (from the TF dataset)
tfs_only <- tf_data$GeneSym

# Find TFs that are present in the gene list (only include TFs that are in the common gene list)
tfs_in_gene_list <- intersect(tfs_only, genes_only)

# Filter interactions where TFs are in the gene list
filtered_tfs_target <- gene_tf_network_data[gene_tf_network_data$target %in% tfs_in_gene_list, ]
filtered_tfs_source <- gene_tf_network_data[gene_tf_network_data$source %in% tfs_in_gene_list, ]

# Merge both filtered interactions
filtered_interactions <- unique(rbind(filtered_tfs_target, filtered_tfs_source))

# Further filter interactions to ensure the connections are only with genes in the common gene list
filtered_interactions <- filtered_interactions[filtered_interactions$target %in% genes_only & filtered_interactions$source %in% genes_only, ]

# Check if there are any interactions
if (nrow(filtered_interactions) == 0) {
  print("No interactions found after filtering. Please check your data.")
} else {
  print(paste("Total interactions after filtering:", nrow(filtered_interactions)))

  # Create the network from the filtered interactions
  gene_tf_graph <- graph_from_data_frame(d = filtered_interactions, directed = TRUE)
  
  # Set colors for TFs and genes (TFs in orange, others in light blue)
  vertex_colors <- ifelse(V(gene_tf_graph)$name %in% tfs_in_gene_list, "orange", "lightblue")

  # Plot the gene-TF interaction network
  plot(gene_tf_graph, vertex.size=5, vertex.label.cex=0.7, vertex.color=vertex_colors, edge.color="gray", main="Filtered Gene-TF Network")
  }



# just for specific Tf_gene network

# Load necessary library
library(igraph)

# Define the highlighted edges (TF-Gene interactions)
highlighted_edges <- data.frame(
  source = c("IKZF1", "IKZF1", "JUN", "JUN", "JUN", "JUN", "JUN", "JUN",
             "RCOR1", "RCOR1", "RCOR1", "RCOR1", "RCOR1", "RCOR1", "RCOR1", "RCOR1",
             "SPI1", "SPI1", "SPI1", "SPI1", "SPI1", "SPI1",
             "ZEB1", "ZEB1", "ZEB1", "ZEB1", "ZEB1", "ZEB1", "ZEB1"),
  target = c("IL15", "PFKFB3", "ACP5", "IL15", "JUN", "PFKFB3", "SPI1", "VEGFA",
             "ACP5", "FLT1", "IL15", "JUN", "PFKFB3", "PFKP", "SPI1", "VEGFA",
             "ACP5", "IL15", "JUN", "PFKP", "SPI1", "VEGFA",
             "ACP5", "IL15", "JUN", "PFKFB3", "PFKP", "SPI1", "VEGFA")
)

# Create graph object
network_graph <- graph_from_data_frame(d = highlighted_edges, directed = TRUE)

# Set node colors: TFs = orange, target genes = light blue
TFs <- unique(highlighted_edges$source)
V(network_graph)$color <- ifelse(V(network_graph)$name %in% TFs, "orange", "lightblue")

# Assign edge width for emphasis
E(network_graph)$width <- 2.5

# Set node size based on degree (more connections = larger size)
V(network_graph)$size <- degree(network_graph, mode = "all") * 2 + 5

# Use a force-directed layout for better spacing
set.seed(123)  # For consistent layout
layout <- layout_with_fr(network_graph)

# Plot the network
plot(
  network_graph,
  layout = layout,
  vertex.label.color = "black",
  vertex.label.cex = 0.8,
  edge.arrow.size = 0.5,
  main = "Transcription Factor - Gene Interaction Network"
)

