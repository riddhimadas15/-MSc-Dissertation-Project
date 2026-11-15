# Limma differential-expression pipeline (RMA normalisation -> limma/eBayes -> filtering -> plots)
# Usage: place the .CEL files for this dataset in ./data_files/ and run from this folder:
#        Rscript limma_GSE93272.R
# Outputs are written to ./output_files/ and ./*.png

#---------------------------------
#	SET PARAMETERS
#---------------------------------
Affy_chip = "hgu133plus2"   # this chip is used for this experiment
gse_current = "GSM2449608"
library(Cairo)
# covdesc file name

#covdesc_file = paste("covdesc",gse_current,sep="_") 
covdesc <- read.delim("covdesc.txt", sep="\t", header=T) 


# Input and output file directories
infile_dir <- "input"
outfile_dir <- "output_files"
dir.create(outfile_dir, showWarnings = FALSE)
#outfile_name <- paste(outfile_dir, "/", gse_current, "_limma_output.csv", sep = "")

test_between <- c("normal_donors", "rheumatoid_arthritis")  # As in covdesc file

#---------------------------------
#	LOAD REQUIRED PACKAGES
#---------------------------------
library(oligo) # Normalization
library(limma) # Differential expression analysis
library(multtest) # Estimating False Discovery Rate (FDR)
library(hgu133plus2.db)
library(pd.hgu133plus2)


#---------------------------------
#	LOAD COVDESC FILE
#---------------------------------
# Read the phenotype description file


colnames(covdesc) <- trimws(colnames(covdesc))
covdesc[] <- lapply(covdesc, trimws)




# Define control and treatment columns
control_cols <- which(covdesc$treatment == test_between[1])
treatment_cols <- which(covdesc$treatment == test_between[2])

control_samples <- length(control_cols)
treatment_samples <- length(treatment_cols)

#---------------------------------
#	LOAD CEL FILES AND NORMALIZE
#---------------------------------
celFiles <- list.celfiles("data_files", full.names = TRUE)
rawData <- read.celfiles(celFiles)

# Add phenotype data to the rawData object
pData(rawData)$gsm <- rownames(pData(rawData))
pData(rawData)$source <- covdesc$treatment

# Boxplot of raw data
xlabs <- as.character(covdesc$source)

# Boxplot for normalized data
Cairo(file = "boxplot.png", 
      type = "png", 
      units = "px", 
      width = 800, 
      height = 600, 
      pointsize = 12, 
      dpi = "auto",
      bg = "white")

# Generate the boxplot
boxplot(rawData, target = "core", boxwex = 0.5, las = 3, names = xlabs)
dev.off()  # Close the device
# Normalize the data using RMA
normalizedData <- rma(rawData)
Cairo(file = "boxplot_2.png", 
      type = "png", 
      units = "px", 
      width = 800, 
      height = 600, 
      pointsize = 12, 
      dpi = "auto",
      bg = "white")
boxplot(normalizedData, target = "core", boxwex = 0.5, las = 3, names = xlabs)
dev.off()
# Get expression values matrix
expval <- exprs(normalizedData)

#---------------------------------
#	DIFFERENTIAL EXPRESSION ANALYSIS USING LIMMA
#---------------------------------
# Create design matrix for the two groups
groups <- pData(rawData)$source
f <- factor(groups, levels = test_between)
design <- model.matrix(~0 + f)
colnames(design) <- test_between

# Fit the model
data.fit <- lmFit(expval, design)

# Create contrasts
contrast.matrix <- makeContrasts(normal_donors - rheumatoid_arthritis, levels = design)

# Fit the contrasts
data.fit.con <- contrasts.fit(data.fit, contrast.matrix)

# Perform moderated t-test
data.fit.eb <- eBayes(data.fit.con)

# Generate volcano plot
Cairo(file = "volcano_plot.png", 
      type = "png", 
      units = "px", 
      width = 800, 
      height = 600, 
      pointsize = 12, 
      dpi = "auto",
      bg = "white")
volcanoplot(data.fit.eb, coef = 1, col = "cornflowerblue")
dev.off()
# Get the top results with FDR correction or Bonferroni 
tab <- topTable(data.fit.eb, coef = 1, number = nrow(normalizedData), adjust.method = "BH", sort.by = "P")
#tab_bonferroni <- topTable(data.fit.eb, coef = 1, number = nrow(normalizedData), adjust.method = "bonferroni", sort.by = "P")

#---------------------------------
#	FOLD CHANGE COMPUTATION
#---------------------------------
# Arithmetic mean fold change
control_data <- rowMeans(2^expval[, control_cols])
treatment_data <- rowMeans(2^expval[, treatment_cols])
FC_linear <- treatment_data / control_data
FC_log2 <- log2(FC_linear)

# Geometric mean fold change
control_geom <- 2^(rowMeans(expval[, control_cols]))
treatment_geom <- 2^(rowMeans(expval[, treatment_cols]))
FC_geom_linear <- treatment_geom / control_geom
FC_geom_log2 <- log2(FC_geom_linear)

# Generate histogram of fold change
Cairo(file = "plot2.png",
      type = "png",
      units = "px",
      width = 600,
      height = 600,
      pointsize = 12,
      dpi = "auto",
      bg = "white")

# Generate the histogram of fold change (Log2)
hist(FC_log2, breaks = 60, xlim = c(-2, 2), main = "Fold Change (Log2)", xlab = "Log2(Fold Change)")

dev.off()  ## Close the device

#---------------------------------
#	ANNOTATION AND OUTPUT
#---------------------------------
annot <- data.frame(
  ENTREZID = sapply(as.list(hgu95av2ENTREZID), paste, collapse = ", "),
  SYMBOL = sapply(as.list(hgu95av2SYMBOL), paste, collapse = ", "),
  DESC = sapply(as.list(hgu95av2GENENAME), paste, collapse = ", ")
)


# Add annotation to normalized data
fData(normalizedData) <- annot[match(rownames(annot), rownames(fData(normalizedData))), ]

# Combine data for output
limma_output <- tab[order(match(rownames(tab), rownames(fData(normalizedData)))), , drop = FALSE]
# Find the intersection of rownames
common_rows <- intersect(rownames(expval), rownames(fData(normalizedData)))
common_rows <- intersect(common_rows, rownames(limma_output))

# Filter each dataset to keep only common rows
expval <- expval[common_rows, , drop = FALSE]
normalizedData <- normalizedData[common_rows, , drop = FALSE]
limma_output <- limma_output[common_rows, , drop = FALSE]

df <- data.frame(
  Gene_ID = rownames(expval),                   # Add gene IDs
  limma_output,                                # Differential expression results
  expval,                                      # Expression values
  Gene_Symbol = fData(normalizedData)$SYMBOL,  # Add gene symbol in a separate column
  Gene_Function = fData(normalizedData)$DESC   # Add gene function in a separate column
)

# Write results to CSV
outfile_files <- "output_files"
#write.csv(df, file = outfile_name, row.names = FALSE)
write.csv(df, paste(outfile_files, "/", gse_current, "limma.csv", sep=""), row.names = FALSE, quote = TRUE)

#---------------------------------
#	FILTERING GENES
#---------------------------------
# Define filtering criteria
# Define filtering criteria
log2FC_cutoff <- 1.0
p_val_cutoff <- 0.05

log2FC_1.5 <- 1.5
p_val_cutoff_1 <- 0.05

# Filter genes based on log2FC and FDR
filtered_genes_1 <- subset(df, 
                         (logFC > log2FC_cutoff | logFC < -log2FC_cutoff) & P.Value < p_val_cutoff)
filtered_genes_1.5 <- subset(df, 
                         (logFC > log2FC_1.5 | logFC < -log2FC_1.5) & P.Value < p_val_cutoff_1)
                                                 
filtered_genes_linear_1.7 <- subset(df, 
                         (logFC > 0.76 | logFC < -0.76) & P.Value < p_val_cutoff_1)
                         
filtered_genes_adj_linear_1.7 <- subset(df, 
                         (logFC > 0.76 | logFC < -0.76) & adj.P.Val < 0.1)

filtered_genes_adj_log_1 <- subset(df, 
                         (logFC > log2FC_cutoff | logFC < -log2FC_cutoff) & adj.P.Val < 0.1)
                         
                         
                         
                                                 
# Subset used for the first heatmap below (|logFC| > 1 and p < 0.05).
# NOTE: the DEG sets reported in the thesis (Table 4.4) use the thresholds listed in the README.
filtered_genes <- filtered_genes_1

# Check how many genes are in the filtered data
message("Number of filtered genes: ", nrow(filtered_genes))

# Check how many genes are in the filtered data
message("Number of filtered genes: ", nrow(filtered_genes_1.5))

# View the filtered genes (first few rows)
print(head(filtered_genes))

# Save filtered genes to a CSV file
filtered_genes_file_1 <- paste(outfile_files, "/", gse_current, "_filtered_genes_1.csv", sep = "")
write.csv(filtered_genes_1, filtered_genes_file_1, row.names = FALSE, quote = TRUE)

filtered_genes_file_1.5 <-paste(outfile_files, "/", gse_current, "_filtered_genes_1.5.csv", sep = "")
write.csv(filtered_genes_1.5, filtered_genes_file_1.5, row.names = FALSE, quote = TRUE)

filtered_genes_file_1.7 <-paste(outfile_files, "/", gse_current, "_filtered_genes_linear_1.7", sep = "")
write.csv(filtered_genes_linear_1.7, filtered_genes_file_1.7, row.names = FALSE, quote = TRUE)

filtered_genes_file_linear_1.7 <-paste(outfile_files, "/", gse_current, "_filtered_genes_adj_linear_1.7", sep = "")
write.csv(filtered_genes_adj_linear_1.7, filtered_genes_file_linear_1.7, row.names = FALSE, quote = TRUE)

filtered_genes_file_adj_log_1 <-paste(outfile_files, "/", gse_current, "_filtered_genes_adj_log_1", sep = "")
write.csv(filtered_genes_adj_log_1, filtered_genes_file_adj_log_1, row.names = FALSE, quote = TRUE)

# Extract expression values for filtered genes
filtered_gene_ids <- rownames(filtered_genes)  # Get gene IDs

# Extract the expression values for filtered genes (ensure the filtered_gene_ids exist in the expval matrix)
filtered_gene_exprs <- expval[filtered_gene_ids, ]

# Ensure numeric data
filtered_gene_exprs <- filtered_gene_exprs <- as.matrix(filtered_gene_exprs)
filtered_gene_exprs <- as.matrix(filtered_gene_exprs)

# Handle NA and Inf values
filtered_gene_exprs <- na.omit(filtered_gene_exprs)  # Remove NA rows
filtered_gene_exprs <- filtered_gene_exprs[apply(filtered_gene_exprs, 1, function(x) all(is.finite(x))), ]

message("Number of filtered_gene_ids: ", length(filtered_gene_ids))
message("Number of matching IDs in expval: ", sum(filtered_gene_ids %in% rownames(expval)))

# Verify if there is any mismatch
missing_ids <- setdiff(filtered_gene_ids, rownames(expval))
if (length(missing_ids) > 0) {
  message("The following gene IDs are missing in expval: ", paste(missing_ids, collapse = ", "))
}
        
# Generate heatmap
library(pheatmap)
library(Cairo)

# Create heatmap plot file

# Generate heatmap
Cairo(file = "filtered_genes_heatmap.png", 
  type = "png", 
  width = 12, 
  height = 10, 
  units = "in", 
  dpi = 300, 
  pointsize = 6)


pheatmap(filtered_gene_exprs, 
     scale = "row", 
     clustering_distance_rows = "correlation", 
     clustering_distance_cols = "correlation", 
     cluster_rows = TRUE, 
     cluster_cols = TRUE, 
     color = colorRampPalette(c("red", "white", "green"))(100), 
     margins = c(8, 8), 
     fontsize = 10, 
     cellwidth = 12, 
     cellheight = 12, 
     show_rownames = TRUE, 
     main = "Heatmap of Filtered Genes")

dev.off()

# Save as PDF
pdf(file = "output_files/filtered_genes_heatmap.pdf", 
    width = 12, height = 20)


pheatmap(filtered_gene_exprs, 
     scale = "row", 
     clustering_distance_rows = "correlation", 
     clustering_distance_cols = "correlation", 
     cluster_rows = TRUE, 
     cluster_cols = TRUE, 
     color = colorRampPalette(c("red", "white", "green"))(100), 
     margins = c(8, 8), 
     fontsize = 10, 
     cellwidth = 12, 
     cellheight = 12, 
     show_rownames = TRUE, 
     main = "Heatmap of Filtered Genes")

dev.off()

message("Heatmap successfully generated.")

#heatmap for filtered genes of logFC 1
# Extract expression values for filtered genes
filtered_gene_ids_1.5 <- rownames(filtered_genes_1.5)  # Get gene IDs

# Extract the expression values for filtered genes (ensure the filtered_gene_ids exist in the expval matrix)
filtered_gene_exprs_1.5 <- expval[filtered_gene_ids_1.5, ]

# Ensure numeric data
filtered_gene_exprs_1.5 <- filtered_gene_exprs_1.5 <- as.matrix(filtered_gene_exprs_1.5)
filtered_gene_exprs_1.5 <- as.matrix(filtered_gene_exprs_1.5)

# Handle NA and Inf values
filtered_gene_exprs_1.5 <- na.omit(filtered_gene_exprs_1.5)  # Remove NA rows
filtered_gene_exprs_1.5 <- filtered_gene_exprs_1.5[apply(filtered_gene_exprs_1.5, 1, function(x) all(is.finite(x))), ]

message("Number of filtered_gene_ids: ", length(filtered_gene_ids_1.5))
message("Number of matching IDs in expval: ", sum(filtered_gene_ids_1.5 %in% rownames(expval)))

# Verify if there is any mismatch
missing_ids_1.5 <- setdiff(filtered_gene_ids_1.5, rownames(expval))
if (length(missing_ids_1.5) > 0) {
  message("The following gene IDs are missing in expval: ", paste(missing_ids_1.5, collapse = ", "))
}
        
# Generate heatmap
library(pheatmap)
library(Cairo)

# Create heatmap plot file

# Generate heatmap
Cairo(file = "filtered_genes_heatmap.png", 
  type = "png", 
  width = 12, 
  height = 10, 
  units = "in", 
  dpi = 300, 
  pointsize = 6)


pheatmap(filtered_gene_exprs_1.5, 
     scale = "row", 
     clustering_distance_rows = "correlation", 
     clustering_distance_cols = "correlation", 
     cluster_rows = TRUE, 
     cluster_cols = TRUE, 
     color = colorRampPalette(c("red", "white", "green"))(100), 
     margins = c(8, 8), 
     fontsize = 10, 
     cellwidth = 12, 
     cellheight = 12, 
     show_rownames = TRUE, 
     main = "Heatmap of Filtered Genes")

dev.off()

# Save as PDF
pdf(file = "output_files/filtered_genes_heatmap_1.5.pdf", 
    width = 12, height = 20)


pheatmap(filtered_gene_exprs_1.5, 
     scale = "row", 
     clustering_distance_rows = "correlation", 
     clustering_distance_cols = "correlation", 
     cluster_rows = TRUE, 
     cluster_cols = TRUE, 
     color = colorRampPalette(c("red", "white", "green"))(100), 
     margins = c(8, 8), 
     fontsize = 10, 
     cellwidth = 12, 
     cellheight = 12, 
     show_rownames = TRUE, 
     main = "Heatmap of Filtered Genes")

dev.off()

message("Heatmap successfully generated.")


# Generate heatmap for unfiltered genes
unfiltered_gene_exprs <- expval  # Or normalizedData
message("Matrix dimensions: ", paste(dim(expval), collapse = " x "))
        
Cairo(file = "unfiltered_genes_heatmap.png", 
      type = "png", 
      width = 15, 
      height = 20, 
      units = "in", 
      dpi = 500, 
      pointsize = 6, 
      compression = "lzw",
      bg = "white")

pheatmap(unfiltered_gene_exprs, 
         scale = "row", 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("red", "white", "green"))(100), 
         margins = c(8, 8), 
         fontsize = 10, 
         cex = 1.0, 
         cellwidth = 12, 
         cellheight = 12, 
         show_rownames = TRUE, 
         main = "Heatmap of Unfiltered Genes")
dev.off()
message("Heatmap successfully generated and saved as 'unfiltered_genes_heatmap.png'.")


