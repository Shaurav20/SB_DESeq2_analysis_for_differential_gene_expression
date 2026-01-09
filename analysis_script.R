## Step 1: Installing and loading required R packages

# Installing required R packages
install.packages(c("pheatmap", "ashr"))

# Checking whether BiocManager has been installed. If not installed, installing it afresh
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Setting working directory
setwd("D:/Projects/New_folder_10")

# Load libraries
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(pheatmap)
library(DEGreport)
library(AnnotationDbi)
library(org.Hs.eg.db)  # or appropriate organism package
library(ashr)  # for apeglm shrinkage
library(RColorBrewer)
library(pheatmap)

## Step 2: Preparing count data

# Read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# Check dimensions
cat("Counts data dimensions:", dim(counts_data), "\n")

# Read in sample info
colData <- read.csv('sample_info.csv')
head(colData)

# Check if the first column of counts_data contains gene identifiers
# If yes, set it as row names
if(colnames(counts_data)[1] == "gene" || colnames(counts_data)[1] == "GeneID") {
  row.names(counts_data) <- counts_data[,1]
  counts_data <- counts_data[,-1]
}

# Making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# Check if they are in the same order
all(colnames(counts_data) == rownames(colData))

# If not in the same order, reorder counts_data columns to match colData
if(!all(colnames(counts_data) == rownames(colData))) {
  cat("Reordering columns to match sample info...\n")
  counts_data <- counts_data[, rownames(colData)]
  all(colnames(counts_data) == rownames(colData))  # Should be TRUE now
}

## Step 3: Construct a DESeqDataSet object

# Check factor levels
cat("Factor levels in dexamethasone:", levels(as.factor(colData$dexamethasone)), "\n")

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)
dds

# Step 3.5: Quality control before filtering

# Plot library sizes
lib_sizes <- colSums(counts(dds))
barplot(lib_sizes, 
        main = "Library sizes", 
        las = 2,
        cex.names = 0.8)

# Plot number of detected genes per sample
detected_genes <- colSums(counts(dds) > 0)
barplot(detected_genes,
        main = "Number of detected genes per sample",
        las = 2,
        cex.names = 0.8)

# Pre-filtering: removing rows with low gene counts
# Improved filtering strategy
cat("Number of genes before filtering:", nrow(dds), "\n")

# Option 1: Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3

# Option 2: More stringent - keep genes with CPM > 1 in at least 2 samples
# library_sizes <- colSums(counts(dds))
# cpm <- counts(dds) / (library_sizes/1e6)
# keep <- rowSums(cpm > 1) >= 2

dds <- dds[keep, ]
cat("Number of genes after filtering:", nrow(dds), "\n")

# Set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

## Step 4: Run DESeq

dds <- DESeq(dds)

## Step 5: Quality Control Plots

# Create a results figures directory
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Plot dispersion estimates
plotDispEsts(dds)

# PCA plot
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
png("results/figures/PCA_plot.png", width = 1200, height = 900)
plotPCA(vsd, intgroup = "dexamethasone") +
  ggtitle("PCA Plot") +
  theme_minimal()
dev.off()

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$dexamethasone, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Ensure clean graphics state
while (!is.null(dev.list())) dev.off()

png("results/figures/Sample_distance_heatmap.png", width = 1200, height = 900)

grid::grid.newpage()

p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows = sampleDists,
              clustering_distance_cols = sampleDists,
              col = colors,
              main = "Sample distance heatmap")

dev.off()

## Step 6: Get Results

res <- results(dds)
res

# Alternative: results with different alpha
res0.01 <- results(dds, alpha = 0.01)
res0.05 <- results(dds, alpha = 0.05)
res0.1  <- results(dds, alpha = 0.1)

# Explore Results ----------------
summary(res)

# Check how many significant genes at different thresholds
cat("\nNumber of significant genes (padj < 0.05):",
    sum(res$padj < 0.05, na.rm = TRUE), "\n")

cat("Number of significant genes (padj < 0.01):",
    sum(res0.01$padj < 0.01, na.rm = TRUE), "\n")

cat("Number of up-regulated genes (padj < 0.05 & log2FC > 1):",
    sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")

cat("Number of down-regulated genes (padj < 0.05 & log2FC < -1):",
    sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")

# Get contrasts
resultsNames(dds)

# Specific contrast
res_treated_vs_untreated <- results(dds,
                                    contrast = c("dexamethasone",
                                                 "treated",
                                                 "untreated"),
                                    alpha = 0.05)

## Step 7: Visualization
           
# MA plot
png("results/figures/MA_plot.png", width = 1200, height = 900)
plotMA(res,
       main = "MA Plot",
       ylim = c(-5, 5))
dev.off()

# Volcano plot

library(ggplot2)

# Prepare data for volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Remove rows with NA adjusted p-values
res_df <- res_df[!is.na(res_df$padj), ]

# Define thresholds
log2FC_cutoff <- 1
padj_cutoff <- 0.05

# Categorize genes
res_df$diffexpressed <- "Not significant"
res_df$diffexpressed[
  res_df$log2FoldChange >= log2FC_cutoff & res_df$padj < padj_cutoff
] <- "Upregulated"

res_df$diffexpressed[
  res_df$log2FoldChange <= -log2FC_cutoff & res_df$padj < padj_cutoff
] <- "Downregulated"

# Plot volcano plot
png("results/figures/Volcano_plot.png", width = 1200, height = 900)
print(
  ggplot(res_df,
         aes(x = log2FoldChange,
             y = -log10(padj))) +
    geom_point(aes(color = diffexpressed),
               alpha = 0.6,
               size = 1.5) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey"
    )) +
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed") +
    labs(title = "Volcano Plot of Differential Expression",
         x = "Log2 fold change",
         y = "-Log10 adjusted p-value") +
    theme_minimal()
)
dev.off()

# Heatmap of top 50 DE genes

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Get significant genes
sig_genes <- res[which(res$padj < 0.05 &
                         abs(res$log2FoldChange) > 1), ]

sig_genes <- sig_genes[order(sig_genes$padj), ]

top_genes <- head(rownames(sig_genes), 50)

# Plot heatmap
if (length(top_genes) > 0) {
  
  heatmap_data <- norm_counts[top_genes, ]
  
  # Ensure clean graphics state
  while (!is.null(dev.list())) dev.off()
  
  png("results/figures/Heatmap_Top50_DEGs.png", width = 1200, height = 1000, res = 150)
  
  # Clear any existing grid graphics
  grid::grid.newpage()
  
  pheatmap(heatmap_data,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           main = "Top 50 Differentially Expressed Genes",
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           annotation_col = colData["dexamethasone"],
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
  dev.off()
}

## Step 8: Gene Annotation

# Add gene symbols if available (example for human)
if (exists("org.Hs.eg.db")) {
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  res_df$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(res_df),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
  
  head(res_df[order(res_df$padj), ], 20)
}

## Step 9: Save Results

# Create results directory
dir.create("results", showWarnings = FALSE)

# Save all results
write.csv(as.data.frame(res), 
          file = "results/DESeq2_results_all.csv",
          row.names = TRUE)

# Save significant results
sig_results <- res[which(res$padj < 0.05), ]
write.csv(as.data.frame(sig_results),
          file = "results/DESeq2_significant.csv",
          row.names = TRUE)

# Save results with gene annotation if available
if(exists("res_df")) {
  write.csv(res_df,
            file = "results/DESeq2_results_annotated.csv",
            row.names = FALSE)
}

# Save normalized counts
write.csv(as.data.frame(norm_counts),
          file = "results/normalized_counts.csv",
          row.names = TRUE)

## Step 10: Generate Report

# Save session information
sink("results/session_info.txt")
sessionInfo()
sink()

# Save summary statistics
sink("results/analysis_summary.txt")
cat("DESeq2 Analysis Summary\n")
cat("=====================\n\n")
cat("Date:", date(), "\n\n")
cat("Input Data:\n")
cat("- Number of samples:", ncol(counts_data), "\n")
cat("- Number of genes before filtering:", nrow(counts_data), "\n")
cat("- Number of genes after filtering:", nrow(dds), "\n\n")
cat("Differential Expression Results:\n")
cat("- Total significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("- Up-regulated (padj < 0.05 & LFC > 1):", sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("- Down-regulated (padj < 0.05 & LFC < -1):", sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")
sink()

cat("\nAnalysis complete! Results saved in 'results/' directory.\n")

