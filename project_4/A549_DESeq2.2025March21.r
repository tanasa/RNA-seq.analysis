library(DESeq2)
library(ggplot2)
library(reshape2)
library(patchwork)  
library(ggrepel)
library(RColorBrewer)
library(glmpca)
library(pheatmap)
library(PoiClaClu)
library(apeglm)
library(ashr)
library(vsn)
library("pheatmap")
library("ReportingTools")

library("BiocParallel")
register(MulticoreParam(4))

# list.files()
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#using-sva-with-deseq2

# Load count data
count_data <- read.csv("A549_featureCounts_table_raw_counts_filter_sum_1000_for_deseq2.txt", 
                       sep="\t", header=TRUE, stringsAsFactors=TRUE)  

head(count_data,2)

rownames(count_data) <- count_data[,1]  # Set Geneid as row names
count_data <- count_data[,-1]           # Remove Geneid column from dataframe

# Print column names to verify correct loading
print(colnames(count_data))  

# Define conditions based on column names
sample_names <- colnames(count_data)

# Assign groups based on sample name prefixes (AC, AP, AV)
conditions <- ifelse(grepl("^AC", sample_names), "AC",
              ifelse(grepl("^AP", sample_names), "AP",
              ifelse(grepl("^AV", sample_names), "AV", NA)))  

# Ensure no unexpected sample names exist
if (any(is.na(conditions))) {
  stop("Some sample names do not match expected patterns (AC, AP, AV). Check column names!")
}

# Create colData ensuring it matches sample names
col_data <- data.frame(
  row.names = sample_names,
  condition = factor(conditions, levels = c("AP", "AC", "AV"))  # Factor levels for correct comparison
)

dim(count_data)
# Remove rows that contain NA values
count_data <- count_data[complete.cases(count_data), ]
dim(count_data)

head(count_data, 2)
tail(count_data, 2)

conditions

# Compute summary statistics for each column
summary_stats <- data.frame(
  Median = apply(count_data, 2, median, na.rm = TRUE),
  Min = apply(count_data, 2, min, na.rm = TRUE),
  Max = apply(count_data, 2, max, na.rm = TRUE)
)

# Print results
print(summary_stats)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "AV")

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)

head(res, 1)
tail(res, 1)

summary(res)
resultsNames(dds)

# rowRanges(dds)
# colData(dds)
# str(dds)

dds <- estimateSizeFactors(dds)
print("the size factors are:")
sizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE) 
head(norm_counts, 2)

write.csv(norm_counts, "A549_normalized_counts.csv", row.names = TRUE)

# Compute summary statistics for each column
summary_stats2 <- data.frame(
  Median = apply(norm_counts, 2, median, na.rm = TRUE),
  Min = apply(norm_counts, 2, min, na.rm = TRUE),
  Max = apply(norm_counts, 2, max, na.rm = TRUE)
)

# Print results
print(summary_stats2, 2)

# perform additional filtering of the normalized counts 

keep <- rowSums(counts(dds, normalized = TRUE) >= 5) >= 4  # Use counts from 'dds'
table(keep)  # Check how many genes will be kept

dds <- dds[keep, ]  
nrow(dds)    # Confirm the number of remaining genes

# Get results for different comparisons
res_AP_vs_AC <- results(dds, contrast = c("condition", "AP", "AC"))
res_AP_vs_AV <- results(dds, contrast = c("condition", "AP", "AV"))
res_AC_vs_AV <- results(dds, contrast = c("condition", "AC", "AV"))

summary(res_AP_vs_AV)
summary(res_AC_vs_AV)
summary(res_AP_vs_AC)

# perform additional filtering of the normalized counts 

# keep <- rowSums(counts(dds, normalized = TRUE) >= 5) >= 4  # Use counts from 'dds'
# table(keep)  # Check how many genes will be kept

# dds <- dds[keep, ]  
# nrow(dds)    # Confirm the number of remaining genes

# Save results
write.csv(as.data.frame(res_AP_vs_AC), file = "A549.DESeq2_AP_vs_AC_results.csv")
write.csv(as.data.frame(res_AP_vs_AV), file = "A549.DESeq2_AP_vs_AV_results.csv")
write.csv(as.data.frame(res_AC_vs_AV), file = "A549.DESeq2_AC_vs_AV_results.csv")

###########################################################
###########################################################
resultsNames(dds)
###########################################################
###########################################################

print("number of differentially bound and expressed transcripts : AP vs AC : pvalue < 0.05")
dim(subset(res_AP_vs_AC, pvalue < 0.05))
dim(subset(res_AP_vs_AC, padj < 0.1))

print("number of differentially bound and expressed transcripts : AP vs AV : pvalue < 0.05")
dim(subset(res_AP_vs_AV, pvalue < 0.05))
dim(subset(res_AP_vs_AV, padj < 0.1))

print("number of differentially bound and expressed transcripts : AC vs AV : pvalue < 0.05")
dim(subset(res_AC_vs_AV, pvalue < 0.05))
dim(subset(res_AC_vs_AV, padj < 0.1))

# type = c("apeglm", "ashr", "normal"),

# If you must use contrast, you should use type="normal" or type="ashr" instead of apeglm, 
# because apeglm only works with coef.  
# Apeglm is the recommended method for log-fold change shrinkage.

# Get results for different comparisons
# resLFCapeglm_AP_vs_AV <- lfcShrink(dds, coef = "condition_AP_vs_AV", type="apeglm")
# resLFCapeglm_AC_vs_AV <- lfcShrink(dds, coef = "condition_AC_vs_AV", type="apeglm")

resLFCashr_AP_vs_AV <- lfcShrink(dds, contrast = c("condition", "AP", "AV"), type="ashr")
resLFCashr_AC_vs_AV <- lfcShrink(dds, contrast = c("condition", "AC", "AV"), type="ashr")
resLFCashr_AP_vs_AC <- lfcShrink(dds, contrast = c("condition", "AP", "AC"), type="ashr")

summary(resLFCashr_AP_vs_AV)
summary(resLFCashr_AC_vs_AV)
summary(resLFCashr_AP_vs_AC)

# Save results
write.csv(as.data.frame(resLFCashr_AP_vs_AC), file = "A549.DESeq2_AP_vs_AC_results.resLFCashr.csv")
write.csv(as.data.frame(resLFCashr_AP_vs_AV), file = "A549.DESeq2_AP_vs_AV_results.resLFCashr.csv")
write.csv(as.data.frame(resLFCashr_AC_vs_AV), file = "A549.DESeq2_AC_vs_AV_results.resLFCashr.csv")

###########################################################
###########################################################
resultsNames(dds)
###########################################################
###########################################################

print("number of differentially bound and expressed transcripts : resLFCashr : AP vs AC : pvalue < 0.05")
dim(subset(resLFCashr_AP_vs_AC, pvalue < 0.05))
dim(subset(resLFCashr_AP_vs_AC, padj < 0.1))

print("number of differentially bound and expressed transcripts : resLFCashr : AP vs AV : pvalue < 0.05")
dim(subset(resLFCashr_AP_vs_AV, pvalue < 0.05))
dim(subset(resLFCashr_AP_vs_AV, padj < 0.1))

print("number of differentially bound and expressed transcripts : resLFCashr : AC vs AV : pvalue < 0.05")
dim(subset(resLFCashr_AC_vs_AV, pvalue < 0.05))
dim(subset(resLFCashr_AC_vs_AV, padj < 0.1))

# Get results for different comparisons

resLFCnormal_AP_vs_AV <- lfcShrink(dds, contrast = c("condition", "AP", "AV"), type="normal")
resLFCnormal_AC_vs_AV <- lfcShrink(dds, contrast = c("condition", "AC", "AV"), type="normal")
resLFCnormal_AP_vs_AC <- lfcShrink(dds, contrast = c("condition", "AP", "AC"), type="normal")

summary(resLFCnormal_AP_vs_AV)
summary(resLFCnormal_AC_vs_AV)
summary(resLFCnormal_AP_vs_AC)

# Save results
write.csv(as.data.frame(resLFCnormal_AP_vs_AC), file = "A549.DESeq2_AP_vs_AC_results.resLFCnormal.csv")
write.csv(as.data.frame(resLFCnormal_AP_vs_AV), file = "A549.DESeq2_AP_vs_AV_results.resLFCnormal.csv")
write.csv(as.data.frame(resLFCnormal_AC_vs_AV), file = "A549.DESeq2_AC_vs_AV_results.resLFCnormal.csv")

###########################################################
###########################################################
resultsNames(dds)
###########################################################
###########################################################

print("number of differentially bound and expressed transcripts : resLFCnormal: AP vs AC : pvalue < 0.05")
dim(subset(resLFCnormal_AP_vs_AC, pvalue < 0.05))
dim(subset(resLFCnormal_AP_vs_AC, padj < 0.1))

print("number of differentially bound and expressed transcripts : resLFCnormal : AP vs AV : pvalue < 0.05")
dim(subset(resLFCnormal_AP_vs_AV, pvalue < 0.05))
dim(subset(resLFCnormal_AP_vs_AV, padj < 0.1))

print("number of differentially bound and expressed transcripts : resLFCashr : AC vs AV : pvalue < 0.05")
dim(subset(resLFCnormal_AC_vs_AV, pvalue < 0.05))
dim(subset(resLFCnormal_AC_vs_AV, padj < 0.1))



# Convert original count data to long format
df_long_before <- melt(as.data.frame(count_data), variable.name = "Sample", value.name = "Expression")
df_long_before$Status <- "Before Normalization" 
# Convert normalized data to long format
df_long_after <- melt(as.data.frame(norm_counts), variable.name = "Sample", value.name = "Expression")
df_long_after$Status <- "After Normalization"  

# Create boxplot for BEFORE normalization
p_before <- ggplot(df_long_before, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Before Normalization", x = "Samples", y = "Expression Levels") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 3000)  

# Create boxplot for AFTER normalization
p_after <- ggplot(df_long_after, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_minimal() +
  labs(title = "After Normalization", x = "Samples", y = "Expression Levels") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 3000)

# Display plots separately (side by side)
p_before + p_after

# Extract raw counts (log-transformed) from DESeq2 object
count_data_log <- log10(count_data + 1)

# Define colors for each sample
colors <- rainbow(ncol(count_data_log))  

# Plot the first sample density
plot(density(count_data_log[, 1]), col = colors[1], lwd = 2,
     main = "before norm : Density Plot for Each Sample",
     xlab = "log10(Counts + 1)", ylab = "Density", ylim = c(0, max(density(count_data_log[, 1])$y)))

# Add density plots for the remaining samples
for (i in 2:ncol(count_data_log)) {
  lines(density(count_data_log[, i]), col = colors[i], lwd = 2)
}

# Add legend
legend("topright", legend = colnames(count_data_log), col = colors, lwd = 2, cex = 0.7)

# Extract normalized counts (log-transformed) from DESeq2 object
count_data2_log <- log10(norm_counts + 1)

# Define colors for each sample
colors <- rainbow(ncol(count_data2_log))  

# Plot the first sample density
plot(density(count_data2_log[, 1]), col = colors[1], lwd = 2,
     main = "after norm : Density Plot for Each Sample",
     xlab = "log10(Counts + 1)", ylab = "Density", 
     ylim = c(0, max(density(count_data_log[, 1])$y)))

# Add density plots for the remaining samples
for (i in 2:ncol(count_data2_log)) {
  lines(density(count_data2_log[, i]), col = colors[i], lwd = 2)
}

# Add legend
legend("topright", legend = colnames(count_data2_log), col = colors, lwd = 2, cex = 0.7)


# Print full summary statistics of the data frame
summary(count_data)

# Compute min, max, and median for each column (sample)
summary_stats <- data.frame(
  Sample = colnames(count_data),
  Min = apply(count_data, 2, min, na.rm = TRUE),
  Max = apply(count_data, 2, max, na.rm = TRUE),
  Median = apply(count_data, 2, median, na.rm = TRUE)
)

# Print the summary table
print(summary_stats)

# Print full summary statistics of the data frame
summary(norm_counts)

# Compute min, max, and median for each column (sample)
summary_stats2 <- data.frame(
  Sample = colnames(norm_counts),
  Min = apply(norm_counts, 2, min, na.rm = TRUE),
  Max = apply(norm_counts, 2, max, na.rm = TRUE),
  Median = apply(norm_counts, 2, median, na.rm = TRUE)
)

# Print the summary table
print(summary_stats2)

dim(norm_counts)
head(norm_counts)

resultsNames(dds)

print("number of differentially bound and expressed transcripts : AP vs AC : pvalue < 0.05")
dim(subset(res_AP_vs_AC, pvalue < 0.05))
dim(subset(res_AP_vs_AC, padj < 0.1))

print("number of differentially bound and expressed transcripts : AP vs AV : pvalue < 0.05")
dim(subset(res_AP_vs_AV, pvalue < 0.05))
dim(subset(res_AP_vs_AV, padj < 0.1))

print("number of differentially bound and expressed transcripts : AC vs AV : pvalue < 0.05")
dim(subset(res_AC_vs_AV, pvalue < 0.05))
dim(subset(res_AC_vs_AV, padj < 0.1))

####################### MA plots

res_AP_vs_AC_filtered <- res_AP_vs_AC[which(res_AP_vs_AC$log2FoldChange > -2 & res_AP_vs_AC$log2FoldChange < 2),]
plotMA(res_AP_vs_AC_filtered, main="AP vs AC", ylim=c(-2,2))

# Filter results where log2FoldChange is between -2 and 2, while removing NAs
resLFCashr_AP_vs_AC_filtered <- na.omit(resLFCashr_AP_vs_AC[which(resLFCashr_AP_vs_AC$log2FoldChange > -2 & 
                                                                 resLFCashr_AP_vs_AC$log2FoldChange < 2), ])
# Plot MA plot for AP vs AC
plotMA(resLFCashr_AP_vs_AC_filtered, main="AP vs AC : ashr shrinkage", ylim=c(-2,2))

####################### MA plots

res_AP_vs_AV_filtered <- res_AP_vs_AV[which(res_AP_vs_AV$log2FoldChange > -2 & res_AP_vs_AV$log2FoldChange < 2),]
plotMA(res_AP_vs_AV, main="AP vs AV", ylim=c(-2,2))

# Corrected filtering for AP vs AV (fixed typo in object name)
resLFCashr_AP_vs_AV_filtered <- na.omit(resLFCashr_AP_vs_AV[which(resLFCashr_AP_vs_AV$log2FoldChange > -2 & 
                                                                 resLFCashr_AP_vs_AV$log2FoldChange < 2), ])
# Plot MA plot for AP vs AV
plotMA(resLFCashr_AP_vs_AV_filtered, main="AP vs AV : ashr shrinkage", ylim=c(-2,2))

####################### MA plots

res_AC_vs_AV_filtered <- res_AC_vs_AV[which(res_AC_vs_AV$log2FoldChange > -2 & res_AC_vs_AV$log2FoldChange < 2),]
plotMA(res_AC_vs_AV, main="AC vs AV", ylim=c(-2,2))

# Corrected filtering for AC vs AV
resLFCashr_AC_vs_AV_filtered <- na.omit(resLFCashr_AC_vs_AV[which(resLFCashr_AC_vs_AV$log2FoldChange > -2 & 
                                                                 resLFCashr_AC_vs_AV$log2FoldChange < 2), ])
# Plot MA plot for AC vs AV
plotMA(resLFCashr_AC_vs_AV_filtered, main="AC vs AV : ashr shrinkage", ylim=c(-2,2))

print("plotting dispersion")

plotDispEsts(dds)
dispersionFunction(dds)

# In order to test for differential expression, we operate on raw counts and use discrete distributions as described 
# in the previous section on differential expression. 
# However for other downstream analyses – e.g. for visualization or clustering – 
# it might be useful to work with transformed versions of the count data.

# Maybe the most obvious choice of transformation is the logarithm. Since count values for a gene can be zero in 
# some conditions (and non-zero in others), some advocate the use of pseudocounts.

# One makes use of the concept of variance stabilizing transformations (VST) (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010),
# and the other is the regularized logarithm or rlog, which incorporates a prior on the sample differences (Love, Huber, and Anders 2014). 
# Both transformations produce transformed data on the log2 scale which has been normalized with respect to library size or 
# other normalization factors.

# Why Are These Transformations Needed?

# RNA-seq count data exhibits a strong mean-variance relationship: genes with low counts tend to have high variance, 
# while genes with high counts tend to have lower relative variance.

# The point of these two transformations, the VST and the rlog, is to remove the dependence of the variance on the mean, 
# particularly the high variance of the logarithm of count data when the mean is low.

# By setting blind to FALSE, the dispersions already estimated will be used to perform transformations, 
# or if not present, they will be estimated using the current design formula.

# Extracting the transformed values : the assay function is used to extract the matrix of normalized values.
# vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)



# Effects of transformations on the variance

rld <- rlog(dds, blind = FALSE)  
vsd <- vst(dds, blind = FALSE) 
ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))
# meanSdPlot(assay(rld))
# meanSdPlot(assay(vsd))

library("pheatmap")

# Select the top 20 differentially expressed genes based on adjusted p-value
top_genes <- rownames(res_AP_vs_AC)[order(res_AP_vs_AC$padj, na.last=NA)][1:20]  #

# Extract normalized transformed counts for the top genes
top_counts <- assay(vsd)[top_genes, ]

# Create annotation dataframe
df <- as.data.frame(colData(dds)["condition"])  # Ensure it is a proper dataframe
colnames(df) <- "Condition"  # Rename column for clarity

# Generate heatmap
pheatmap(top_counts, 
         cluster_rows=TRUE,  # Cluster rows to group similar genes
         show_rownames=TRUE,  # Show gene names
         cluster_cols=TRUE,  # Cluster samples
         annotation_col=df,  # Add sample condition annotations
         scale="row",  # Normalize each gene's expression across samples
         fontsize_row=8)  # Adjust row text size for readability



print("PCA plot of rld-transformed counts")

rld <- rlog(dds, blind = FALSE)  # Regularized log transformation
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  scale_color_brewer(palette="Dark2") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()

print("MDS plot of rld-transformed counts")

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
mdsData <- cmdscale(sampleDistMatrix)
mds_df <- data.frame(MDS1 = mdsData[,1], MDS2 = mdsData[,2], condition = col_data$condition)

ggplot(mds_df, aes(MDS1, MDS2, color=condition)) +
  geom_point(size=4) +
  scale_color_brewer(palette="Dark2") +
  ggtitle("MDS Plot of Samples") +
  theme_minimal()

print("PCA plot of vst-transformed counts")

vsd <- vst(dds, blind = FALSE) 
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  scale_color_brewer(palette="Dark2") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()

# MDS Plot
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
mdsData <- cmdscale(sampleDistMatrix)
mds_df <- data.frame(MDS1 = mdsData[,1], MDS2 = mdsData[,2], condition = col_data$condition)

ggplot(mds_df, aes(MDS1, MDS2, color=condition)) +
  geom_point(size=4) +
  scale_color_brewer(palette="Dark2") +
  ggtitle("MDS Plot of Samples") +
  theme_minimal()

# str(dds)

print("PCA by using GLMPCA library. RLOG and VSD transformation are more suitable than SCALE() function.")
options(repr.plot.width = 10, repr.plot.height = 14)

gpca <- glmpca(assay(dds), L=2)
gpca.dat <- gpca$factors

# Ensure gpca.dat has correct sample metadata
gpca.dat$sample <- colnames(dds)              # Assign sample names if missing
gpca.dat$condition <- colData(dds)$condition  # Assign condition/group labels

str(gpca.dat)
head(gpca.dat, 2)

# **Improved GLM-PCA Plot with Sample Labels**
p <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition)) +
  geom_point(size = 6, alpha = 0.8, stroke = 1) +  
  geom_text_repel(aes(label = sample), size = 5, box.padding = 0.5, max.overlaps = 10) +  
  coord_fixed() +  
  theme_minimal(base_size = 18) +  
  labs(title = "GLM-PCA on Filtered Counts", 
       x = "GLM-PC1", y = "GLM-PC2") +
  scale_color_brewer(palette = "Set2") +  
  theme(legend.position = "right",
        plot.title = element_text(size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) +
  guides(shape = guide_legend(override.aes = list(size = 5)))

p

print("perform GLM-PCA")

# Extract count matrix from dds
dds_matrix <- assay(dds)
dim(dds_matrix)

# Remove rows (genes) that are all zeros
dds_matrix_filtered <- dds_matrix[rowSums(dds_matrix) > 0, ]
dim(dds_matrix_filtered)

# Run GLM-PCA on the filtered count data
gpca <- glmpca(dds_matrix_filtered, L = 2)
gpca.dat <- as.data.frame(gpca$factors)
gpca.dat$sample <- colnames(dds)

# Assign condition labels from metadata
if ("condition" %in% colnames(colData(dds))) {
  gpca.dat$condition <- colData(dds)$condition
} else {
  stop("Column 'condition' not found in colData(dds). Check available columns.")
}

p <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition)) +
  geom_point(size = 6, alpha = 0.8, stroke = 1) +  
  geom_text_repel(aes(label = sample), size = 5, box.padding = 0.5, max.overlaps = 10) +  
  coord_fixed() +  
  theme_minimal(base_size = 18) +  
  labs(title = "GLM-PCA on raw, not-norm counts", 
       x = "GLM-PC1", y = "GLM-PC2") +
  scale_color_brewer(palette = "Set2") +  
  theme(legend.position = "right",
        plot.title = element_text(size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) +
  guides(shape = guide_legend(override.aes = list(size = 5)))  

print(p)

print("Heatmap of RLOG transformed matrix") 

rlog_matrix <- assay(rld)

# Compute Euclidean distance
sampleDists <- dist(t(rlog_matrix))
sampleDistMatrix <- as.matrix(sampleDists)

# Define improved color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate Heatmap**
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         fontsize_row = 10,
         fontsize_col = 10,
         cellwidth = 60,
         cellheight = 60,
         angle_col = 45,
         main = "Sample Distance Heatmap (RLOG)")

print("Heatmap of VST transformed matrix") 

vsd_matrix <- assay(vsd)  

# Compute Euclidean distance
sampleDists <- dist(t(vsd_matrix))  
sampleDistMatrix <- as.matrix(sampleDists)  

# Define an improved color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         fontsize_row = 10,  # Make labels readable
         fontsize_col = 10,
         cellwidth = 60,  # Adjust cell size
         cellheight = 60,
         angle_col = 45,  # Rotate column labels for readability
         main = "Sample Distance Heatmap (VST)")

# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package.
# This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating
# the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows
# instead of columns, so we need to transpose the counts in dds.

library("PoiClaClu")

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)

# Extract correct sample names
sample_names <- rownames(colData(dds))  
if (is.null(sample_names)) {
  sample_names <- colnames(dds)
}

# Assign row and column names to keep them in the heatmap
rownames(samplePoisDistMatrix) <- sample_names
colnames(samplePoisDistMatrix) <- sample_names

# Define an improved color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# **Generate Heatmap with Proper Labels**
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         fontsize_row = 10,  # Make labels readable
         fontsize_col = 10,
         cellwidth = 60,  # Adjust cell size
         cellheight = 60,
         angle_col = 45,  # Rotate column labels for readability
         main = "Poisson Distance Heatmap on raw, not-norm counts")

# str(dds)

# library("ReportingTools")
# res = results(dds, contrast=c("condition", "AP", "AC"))
# tmp <- tempdir(getwd()) # you would instead use a meaningful path here
# rep <- HTMLReport(shortName = paste("AP.vs.AC", sep="."),  
#                  title = ,
#                  basePath = tmp,
#                  reportDirectory = "report")
# publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
# finish(rep)

# library("ReportingTools")
# res = results(dds, contrast=c("condition", "AP", "AV"))
# tmp <- tempdir(getwd()) # you would instead use a meaningful path here
# rep <- HTMLReport(shortName = paste("AP.vs.AV", sep="."),  
#                  title = ,
#                  basePath = tmp,
#                  reportDirectory = "report")
# publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
# finish(rep)

# library("ReportingTools")
# res = results(dds, contrast=c("condition", "AC", "AV"))
# tmp <- tempdir(getwd()) # you would instead use a meaningful path here
# rep <- HTMLReport(shortName = paste("AC.vs.AV", sep="."),  
#                  title = ,
#                  basePath = tmp,
#                  reportDirectory = "report")
# publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
# finish(rep)

print("SVA analysis")
# SV1, SV2, ... are surrogate variables — latent (hidden) factors estimated from the data that capture unwanted variation 
# (like batch effects, technical noise, or hidden biological subtypes).
# You can think of them as "virtual covariates" — constructed purely from the structure of your data — 
# that explain sources of variation not included in your model (like treatment or condition).

library(sva)
library(RUVSeq)

############################################################################################################################### using SVAseq

dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv = 2)
head(svseq$sv, 2)

par(mfrow = c(2, 1),           # 2 rows, 1 column layout
    mar = c(4, 5, 4, 2) + 0.1, # bottom, left, top, right margins
    cex.main = 1.2,            # title size
    cex.axis = 1,              # axis label size
    cex.lab = 1.1)             # axis title size

for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$condition,
             vertical = TRUE,
             method = "jitter",
             pch = 16,
             col = "steelblue",
             main = paste0("Surrogate Variable SV", i),
             ylab = "SV Value",
             xlab = "Condition")
  abline(h = 0, lty = 2, col = "gray")
}


# Finally, in order to use SVA to remove any effect on the counts from our surrogate variables, we simply add these two surrogate variables 
# as columns to the DESeqDataSet and then add them to the design:

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + condition
  
ddssva$SV1
ddssva$SV2

# length(ddssva$SV1)
# length(ddssva$SV2)

ddssva <- DESeq(ddssva)
resultsNames(ddssva)

# rowRanges(ddssva)
# colData(ddssva)
# assays(ddssva)
# assay(ddssva)
# length(rowRanges(ddssva))

res_ddssva <- results(ddssva)
resultsNames(res_ddssva)

# Get results for different comparisons
res_ddssva_AP_vs_AC <- results(ddssva, contrast = c("condition", "AP", "AC"))
res_ddssva_AP_vs_AV <- results(ddssva, contrast = c("condition", "AP", "AV"))
res_ddssva_AC_vs_AV <- results(ddssva, contrast = c("condition", "AC", "AV"))

summary(res_ddssva_AP_vs_AV)
summary(res_ddssva_AC_vs_AV)
summary(res_ddssva_AP_vs_AC)

# Save results
write.csv(as.data.frame(res_ddssva_AP_vs_AC), file = "A549.DESeq2_AP_vs_AC_results.sva.csv")
write.csv(as.data.frame(res_ddssva_AP_vs_AV), file = "A549.DESeq2_AP_vs_AV_results.sva.csv")
write.csv(as.data.frame(res_ddssva_AC_vs_AV), file = "A549.DESeq2_AC_vs_AV_results.sva.csv")

###########################################################
###########################################################

print("number of differentially bound and expressed transcripts : AP vs AC : pvalue < 0.05")
dim(subset(res_ddssva_AP_vs_AC, pvalue < 0.05))
dim(subset(res_ddssva_AP_vs_AC, padj < 0.1))

print("number of differentially bound and expressed transcripts : AP vs AV : pvalue < 0.05")
dim(subset(res_ddssva_AP_vs_AV, pvalue < 0.05))
dim(subset(res_ddssva_AP_vs_AV, padj < 0.1))

print("number of differentially bound and expressed transcripts : AC vs AV : pvalue < 0.05")
dim(subset(res_ddssva_AC_vs_AV, pvalue < 0.05))
dim(subset(res_ddssva_AC_vs_AV, padj < 0.1))

########################################################################################################################
########################################################################################################################

options(repr.plot.width = 6, repr.plot.height = 6)
vsd2 <- vst(ddssva, blind=TRUE)
plotPCA(vsd2, "condition")

options(repr.plot.width = 6, repr.plot.height = 6)
rld2 <- rlog(ddssva, blind = TRUE)
plotPCA(rld2, "condition")

# We would then run DESeq with the new design to re-estimate the parameters and results.

print("RUVseq analysis")

set <- newSeqExpressionSet(counts(dds))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(res)[which(res$pvalue > .1)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)
pData(set)

par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$condition, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
 }

ddsruv <- dds
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + condition

# We would then run DESeq with the new design to re-estimate the parameters and results.
# length(ddsruv$SV1)
# length(ddsruv$SV2)

ddsruv <- DESeq(ddsruv)
resultsNames(ddsruv)

# rowRanges(ddsruv)
# colData(ddsruv)
# assays(ddsruv)
# assay(ddsruv)
# length(rowRanges(ddsruv))

res_ddsruv <- results(ddsruv)
resultsNames(res_ddsruv)

# Get results for different comparisons
res_ddsruv_AP_vs_AC <- results(ddsruv, contrast = c("condition", "AP", "AC"))
res_ddsruv_AP_vs_AV <- results(ddsruv, contrast = c("condition", "AP", "AV"))
res_ddsruv_AC_vs_AV <- results(ddsruv, contrast = c("condition", "AC", "AV"))

summary(res_ddsruv_AP_vs_AV)
summary(res_ddsruv_AC_vs_AV)
summary(res_ddsruv_AP_vs_AC)

# Save results
write.csv(as.data.frame(res_ddsruv_AP_vs_AC), file = "A549.DESeq2_AP_vs_AC_results.ruv.csv")
write.csv(as.data.frame(res_ddsruv_AP_vs_AV), file = "A549.DESeq2_AP_vs_AV_results.ruv.csv")
write.csv(as.data.frame(res_ddsruv_AC_vs_AV), file = "A549.DESeq2_AC_vs_AV_results.ruv.csv")

###########################################################
###########################################################

print("number of differentially bound and expressed transcripts : AP vs AC : pvalue < 0.05")
dim(subset(res_ddsruv_AP_vs_AC, pvalue < 0.05))
dim(subset(res_ddsruv_AP_vs_AC, padj < 0.1))

print("number of differentially bound and expressed transcripts : AP vs AV : pvalue < 0.05")
dim(subset(res_ddsruv_AP_vs_AV, pvalue < 0.05))
dim(subset(res_ddsruv_AP_vs_AV, padj < 0.1))

print("number of differentially bound and expressed transcripts : AC vs AV : pvalue < 0.05")
dim(subset(res_ddsruv_AC_vs_AV, pvalue < 0.05))
dim(subset(res_ddsruv_AC_vs_AV, padj < 0.1))


