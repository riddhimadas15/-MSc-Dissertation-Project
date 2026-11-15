# Load required libraries
library(pheatmap)
library(Cairo)  # For high-quality PNG output

# Read the data
data <- read.csv("common.csv", header = TRUE)

# Extract expression values
expression_data <- data[, sapply(data, is.numeric)]
# Convert to matrix
expression_matrix <- as.matrix(expression_data)

# Set row names as Gene Symbols instead of Gene_ID
rownames(expression_matrix) <- data$SYMBOL

# Save the heatmap as a high-quality PNG with improved readability
CairoPNG("expression_plot_66.png", width = 3000, height = 2500, res = 300)

# Generate heatmap with Gene Symbols as row names
pheatmap(expression_matrix, 
         scale = "row",  # Standardize rows (genes) to show relative expression
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE,  # Display Gene Symbols
         fontsize_row = 6)  # Reduce font size for readability

# Close the device
dev.off()

# Print confirmation message
cat("Heatmap saved as 'expression_plot_66.png' with Gene Symbols\n")
