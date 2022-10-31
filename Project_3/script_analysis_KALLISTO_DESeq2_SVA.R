################################################################################################################################ 
################################################################################################################################
################################################################################################################################
library("ggplot2")
library("reshape2")
library("data.table")
library("DESeq2")
library("limma")
library("Glimma")
library("edgeR")
library("pheatmap")
library("ComplexHeatmap")
library("scatterplot3d")
library("enrichR")
library("tidyr")
library("plyr")
library("dplyr")
library("RColorBrewer")
# library("VennDiagram") 
# library("Vennerable")
library("gplots")
library("scatterplot3d")
library("cowplot")
library("tximport")
library("tximeta")
library("pheatmap")
library("RColorBrewer")
library("ReportingTools")
library("Glimma")
library("sva")
library("RUVSeq")
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
######################################################################################################### NAME="analysis.DEseq2"

NAME = "KALLISTO.RefSeq.INTEGRATED.TPM.txt.no.dup.genes.txt"

counts = read.delim("KALLISTO.RefSeq.INTEGRATED.TPM.txt.no.dup.genes.txt", 
                     sep="\t", header=T, stringsAsFactors=F)

rownames(counts) = counts$gene_name

colnames(counts)
# [1] "target_id.S1"  "length.S1"     "eff_length.S1" "est_counts.S1"
# [5] "tpm.S1"        "SYMBOL.x"      "length.S2"     "eff_length.S2"
# [9] "est_counts.S2" "tpm.S2"        "SYMBOL.y"      "length.S3"    
#[13] "eff_length.S3" "est_counts.S3" "tpm.S3"        "SYMBOL.x.1"   
#[17] "length.S4"     "eff_length.S4" "est_counts.S4" "tpm.S4"       
#[21] "SYMBOL.y.1"    "length.S5"     "eff_length.S5" "est_counts.S5"
#[25] "tpm.S5"        "SYMBOL.x.2"    "length.S6"     "eff_length.S6"
#[29] "est_counts.S6" "tpm.S6"        "SYMBOL.y.2"    "length.S7"    
#[33] "eff_length.S7" "est_counts.S7" "tpm.S7"        "SYMBOL.x.3"   
#[37] "length.S8"     "eff_length.S8" "est_counts.S8" "tpm.S8"       
#[41] "SYMBOL.y.3"    "length.S9"     "eff_length.S9" "est_counts.S9"
#[45] "tpm.S9"        "SYMBOL"

samples = read.delim("KALLISTO.RefSeq.INTEGRATED.TPM.txt.no.dup.genes.txt.COUNTS.txt.the.SAMPLE.NAMES.txt", 
                     sep="\t", header=T, stringsAsFactors=F)

samples

#  sample  group     experiment batch
#1     S1 lacCTR   lacZshRNA-WT     1
#2     S2 lacCTR   lacZshRNA-WT     2
#3     S3 lacSTZ  lacZshRNA-STZ     3
#4     S4 lacSTZ  lacZshRNA-STZ     2
#5     S5 zfpCTR  Zfp36shRNA-WT     4
#6     S6 zfpCTR  Zfp36shRNA-WT     4
#7     S7 zfpSTZ Zfp36shRNA-STZ     5
#8     S8 zfpSTZ Zfp36shRNA-STZ     5
#9     S9 zfpSTZ Zfp36shRNA-STZ     5

# in this version of the analysis : 

samples$sample = as.factor(samples$sample)
samples$batch = as.factor(samples$batch)
samples$group = as.factor(samples$group)
samples$experiment = as.factor(samples$experiment)

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ TPM
################################################################################################################################ 
################################################################################################################################ 

counts.TPM <- subset(counts, select=c("SYMBOL", 
                "tpm.S1",
                "tpm.S2",
                "tpm.S3", 
                "tpm.S4", 
                "tpm.S5",
                "tpm.S6", 
                "tpm.S7", 
                "tpm.S8"))

head(counts.TPM)
dim(counts.TPM)  
rownames(counts.TPM) = counts.TPM$SYMBOL

################################################################################################################################ 
################################################################################################################################

counts.SMALL.TPM = counts.TPM 
dim(counts.SMALL.TPM)

################################################################################################################################ 
################################################################################################################################

counts.TPM = counts.TPM[,-1]

### changing the names of the columns :

colnames(counts.TPM)[1] = "S1"
colnames(counts.TPM)[2] = "S2"
colnames(counts.TPM)[3] = "S3"
colnames(counts.TPM)[4] = "S4"
colnames(counts.TPM)[5] = "S5"
colnames(counts.TPM)[6] = "S6"
colnames(counts.TPM)[7] = "S7"
colnames(counts.TPM)[8] = "S8"
# colnames(counts.TPM)[9] = "S9"

################################################################################################################################
################################################################################################################################
# colnames(counts.TPM)
# pch <- c(0, 1, 2, 5, 15, 16, 17, 18, 19)
# [1] "S1" "S2" "S3" "S4" "S5" "S6" "S7" "S8" "S9"
# colors <- c("darkgreen", "darkgreen", "red", "red", "brown", "brown", "blue", "blue", "blue")
################################################################################################################################
################################################################################################################################
############################################################################################################# displaying the PCA

pca <- prcomp(t(counts.TPM))

### plotting the PCA :

pdf(paste(NAME, "PCA", "on.TPM", "in.3D.pdf", sep="."), width=10, height=10)

s3d <- scatterplot3d(pca$x[,1:3],
                           color = c(rep("darkgreen",2), rep("red",2), rep("brown",2), rep("blue",2)),
                           # pch=18,
                           # pch = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
                           pch = c(0, 1, 2, 5, 15, 16, 17, 18),
                           main = "PCA analysis (TPM)", 
                           grid = TRUE, 
                           box = TRUE)

legend("topright",
        s3d$xyz.convert(18, 0, 12), 
        # s3D$xyz.convert(7.5, 3, 4.5), 
        legend = row.names(pca$x[,1:3]),
        col = c(rep("darkgreen",2), rep("red",2), rep("brown",2), rep("blue",2)),
        # pch=18,
        # pch = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
        pch = c(0, 1, 2, 5, 15, 16, 17, 18))

dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ RAW
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
### another WAY to change the NAMES of COLUMNS

counts.small <- subset(counts, select=c("SYMBOL", 
                "est_counts.S1",
                "est_counts.S2",
                "est_counts.S3", 
                "est_counts.S4", 
                "est_counts.S5",
                "est_counts.S6", 
                "est_counts.S7", 
                "est_counts.S8"))

head(counts.small)
dim(counts.small) 

############################################################################################# the 1st COLUMN is the GENE SYMBOL 
########################################################################################### changing the names of the columns :

colnames(counts.small)[2] = "S1"
colnames(counts.small)[3] = "S2"
colnames(counts.small)[4] = "S3"
colnames(counts.small)[5] = "S4"
colnames(counts.small)[6] = "S5"
colnames(counts.small)[7] = "S6"
colnames(counts.small)[8] = "S7"
colnames(counts.small)[9] = "S8"

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

write.table(as.data.frame(counts.small),
            file=paste(NAME, "COUNTS.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

### the ROW NAMES :

rownames(counts.small) = counts.small$SYMBOL
counts.small = counts.small[,-1]

#  head(samples)
#  sample  group    experiment batch
#1     S1 lacCTR  lacZshRNA-WT     1
#2     S2 lacCTR  lacZshRNA-WT     2
#3     S3 lacSTZ lacZshRNA-STZ     3
#4     S4 lacSTZ lacZshRNA-STZ     2
#5     S5 zfpCTR Zfp36shRNA-WT     4
#6     S6 zfpCTR Zfp36shRNA-WT     4

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
############################################################################################################# displaying the PCA

pca <- prcomp(t(counts.small))

### plotting the PCA :

pdf(paste(NAME, "PCA", "on.RAW", "in.3D.pdf", sep="."), width=10, height=10)

s3d <- scatterplot3d(pca$x[,1:3],
                           color = c(rep("darkgreen",2), rep("red",2), rep("brown",2), rep("blue",2)),
                           # pch=18,
                           # pch = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
                           pch = c(0, 1, 2, 5, 15, 16, 17, 18),
                           main = "PCA analysis (RAW)", 
                           grid = TRUE, 
                           box = TRUE)

legend("topright",
        s3d$xyz.convert(18, 0, 12), 
        # s3D$xyz.convert(7.5, 3, 4.5), 
        legend = row.names(pca$x[,1:3]),
        col = c(rep("darkgreen",2), rep("red",2), rep("brown",2), rep("blue",2)),
        # pch=18,
        # pch = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
        pch = c(0, 1, 2, 5, 15, 16, 17, 18))

dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ DESeq2 
################################################################################################################################
########################################################################################### starting to use the CODE from DESeq2 

head(counts.small)
samples

################################################################################################################################ 
################################################################################################################################ SUMMARIZED EXPERIMENT 
################################################################################################################################

gse = SummarizedExperiment(assays=list(counts=as.matrix(counts.small)), colData=samples)

gse
assayNames(gse)
rowRanges(gse)
colData(gse)
colSums(assay(gse))

gse$sample
gse$group
gse$experiment
gse$batch

# gse$sample
#[1] "S1" "S2" "S3" "S4" "S5" "S6" "S7" "S8" "S9"
# gse$group
#[1] "lacCTR" "lacCTR" "lacSTZ" "lacSTZ" "zfpCTR" "zfpCTR" "zfpSTZ" "zfpSTZ"
#[9] "zfpSTZ"
# gse$experiment
#[1] "lacZshRNA-WT"   "lacZshRNA-WT"   "lacZshRNA-STZ"  "lacZshRNA-STZ" 
#[5] "Zfp36shRNA-WT"  "Zfp36shRNA-WT"  "Zfp36shRNA-STZ" "Zfp36shRNA-STZ"
#[9] "Zfp36shRNA-STZ"
# gse$batch

# gse$sample = as.factor(gse$sample)
# gse$batch = as.factor(gse$batch)
# gse = SummarizedExperiment(assays=list(counts=as.matrix(counts.small)), colData=samples)

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################

# In DESeq2, the custom class is called DESeqDataSet. It is built on top of the SummarizedExperiment class, 
# and it is easy to convert SummarizedExperiment objects into DESeqDataSet objects, which we show below. 
# One of the two main differences is that the assay slot is instead accessed using the counts accessor 
# function, and the DESeqDataSet class enforces that the values in this matrix are non-negative integers.

# A second difference is that the DESeqDataSet has an associated design formula. 
# The experimental design is specified at the beginning of the analysis, 
# as it will inform many of the DESeq2 functions how to treat the samples in the analysis.

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 

dds <- DESeqDataSetFromMatrix( countData = round(as.matrix(counts.small)),
                               colData=samples, 
                               design = ~ group)

###

rowRanges(dds)
colData(dds)
assays(dds)
assay(dds)

length(rowRanges(dds)

png(paste(NAME, "part1.plot.DENSITIES.before.FILTERING.png", sep="."))
boxplot(log10(counts(dds)+1), legend = TRUE, main = "before filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ FILTERING

# a minimal filtering to reduce the size of the dataset. We do not need to retain genes if they do not have a 
# count of 5 or more for 4 or more samples as these genes will have no statistical power to detect differences

keep <- rowSums(counts(dds) >= 5) >= 4
table(keep)
#  keep
#  FALSE  TRUE 
#  4503 12927

dds <- dds[keep,]

length(rowRanges(dds))
# [1] 12927

nrow(dds)

boxplot(log10(counts(dds)+1))

# Visualize the distribution of gene expression levels AFTER FILTERING

png(paste(NAME, "part1.plot.DENSITIES.after.FILTERING.png", sep="."))
boxplot(log10(counts(dds)+1), legend = TRUE, main = "after filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ NORMALIZATION

# The main function in DESeq2 involves the computation of the size factors which normalize for differences in sequencing depth among samples. 
# We can also compute these size factors manually, so that the normalized counts are available for plotting:

dds <- estimateSizeFactors(dds)

png(paste(NAME, "part1.plot.DENSITIES.after.NORMALIZATION.png", sep="."))
boxplot(log10(counts(dds, normalized=TRUE) + 1))
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ TRANSFORMATION
################################################################################################################################ 
################################################################################################################################ EDA and VST
################################################################################################################################ 

# Taking the logarithm of counts plus a pseudocount of 1 is a common transformation, 
# but it tends to inflate the sampling variance of low counts such that it is even larger than biological variation across groups of samples. 
# In DESeq2 we therefore provide transformations which produce log-scale data such that the systematic trends have been removed. 
# Our recommended transformation is the variance-stabilizing transformation, or VST, and it can be called with the vst function:

vsd <- vst(dds, blind=TRUE)

class(vsd)
assay(vsd)[1:3,1:3]
colData(vsd)

# The VST data is appropriate for calculating distances between samples or for performing PCA

png(paste(NAME, "part2.VST.PCA.group.png", sep="."))
plotPCA(vsd, "group")
dev.off()

png(paste(NAME, "part2.VST.PCA.sample.png", sep="."))
plotPCA(vsd, "sample")
dev.off()

png(paste(NAME, "part2.VST.PCA.group.samples.v1.png", sep="."))
pcaData <- plotPCA(vsd, intgroup = c("group", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape = sample)) +
  geom_point(shape = c(1,2,3,4,5,6,7,8), size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text(
    label = samples$sample, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )
dev.off()

png(paste(NAME, "part2.VST.PCA.group.samples.v2.png", sep="."))
pcaData <- plotPCA(vsd, intgroup = c("group", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = group, shape = sample)) +
  geom_point(shape = c(1,2,3,4,5,6,7,8), size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_label(
    label = samples$sample, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )
dev.off()

################################################################################################################################
################################################################################################################################ 
################################################################################################################################

# Note that we do not recommend working with the transformed data for the primary differential expression analysis. 
# Instead we will use the original counts and a generalized linear model (GLM) which takes into account the expected variance from either low or high counts. 
# For statistical details, please refer to the DESeq2 methods paper (Love, Huber, and Anders 2014)

################################################################################################################################ TRANSFORMATION
################################################################################################################################ 
################################################################################################################################ EDA and RLOG
################################################################################################################################ 
######################################################## For a fully unsupervised transformation, one can set blind = TRUE (which is the default) 

rld <- rlog(dds, blind = TRUE)
head(assay(rld), 3)
colData(rld)

# The RLD data : 

png(paste(NAME, "part2.RLD.PCA.group.png", sep="."))
plotPCA(rld, "group")
dev.off()

png(paste(NAME, "part2.RLD.PCA.sample.png", sep="."))
plotPCA(rld, "sample")
dev.off()

png(paste(NAME, "part2.RLD.PCA.group.samples.v1.png", sep="."))
pcaData <- plotPCA(rld, intgroup = c("group", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = sample)) +
  geom_point(shape = c(1,2,3,4,5,6,7,8), size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  geom_text(
    label = samples$sample, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T)
dev.off()

png(paste(NAME, "part2.RLD.PCA.group.samples.v2.png", sep="."))
pcaData <- plotPCA(rld, intgroup = c("group", "sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = sample)) +
  geom_point(shape = c(1,2,3,4,5,6,7,8), size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  geom_label(
    label = samples$sample, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T)
dev.off()

################################################################################################################################ TRANSFORMATION
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ using GLM-PCA
################################################################################################################################

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors

gpca.dat$group <- dds$group
gpca.dat$sample <- dds$sample

png(paste(NAME, "part2.using.GLM.PCA.group.samples.v1.png", sep="."))
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, shape = sample)) +
                 # geom_point(size =3) +  
                 geom_point(shape = c(1,2,3,4,5,6,7,8), size = 5) +
                 coord_fixed() + 
                 ggtitle("glmpca - Generalized PCA") + 
                 geom_text(
                         label = samples$sample, 
                         nudge_x = 0.25, nudge_y = 0.25, 
                 check_overlap = T)
dev.off()

png(paste(NAME, "part2.using.GLM.PCA.group.samples.v2.png", sep="."))
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group, shape = sample)) +
                 # geom_point(size =3) +  
                 geom_point(shape = c(1,2,3,4,5,6,7,8), size = 5) +
                 coord_fixed() + 
                 ggtitle("glmpca - Generalized PCA") + 
                 geom_label(
                        label = samples$sample, 
                        nudge_x = 0.25, nudge_y = 0.25, 
                 check_overlap = T)
dev.off()

################################################################################################################################ 
################################################################################################################################ SAMPLE DISTANCES
################################################################################################################################ 
# We use the R function dist to calculate the Euclidean distance between samples. To ensure we have a roughly equal contribution from all genes, 
# we use it on the VST data.
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, 
# and different dimensions (here, genes) to be columns.

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(paste(NAME, "part3.SAMPLE.DISTANCES.VSD.euclidean.png", sep="."))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. 
# This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating 
# the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows 
# instead of columns, so we need to transpose the counts in dds.

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix(poisd$dd )
rownames(samplePoisDistMatrix) <- paste(dds$sample, dds$group, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL

png(paste(NAME, "part3.SAMPLE.DISTANCES.VSD.poisson.png", sep="."))
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()

################################################################################################################################ 
################################################################################################################################ MDS plot using VST
################################################################################################################################ 
################################################################################################################################

mds <- as.data.frame(colData(vsd))  %>%  cbind(cmdscale(sampleDistMatrix))

png(paste(NAME, "part4.VST.MDS.euclidean.png", sep="."), width = 480, height = 480)
ggplot(mds, aes(x = `1`, y = `2`, color = group, shape = sample )) +
            # geom_point(size = 3) + 
            geom_point(shape = c(1,2,3,4,5,6,7,8), size = 8) +
            coord_fixed() + 
            ggtitle("MDS with VST data : euclidean") + 
            geom_text(    
                         label = samples$sample, 
                         nudge_x = 0.25, nudge_y = 0.25, 
                         check_overlap = T)
dev.off()

mdsPois <- as.data.frame(colData(dds)) %>% cbind(cmdscale(samplePoisDistMatrix))

png(paste(NAME, "part4.VST.MDS.poisson.png", sep="."), width = 480, height = 480)
ggplot(mdsPois, aes(x = `1`, y = `2`, color = group, shape = sample)) +
            # geom_point(size = 3) +   
            geom_point(shape = c(1,2,3,4,5,6,7,8), size = 8) +
            coord_fixed() + 
            ggtitle("MDS with VST data : poisson") + 
            geom_text(    
                         label = samples$sample, 
                         nudge_x = 0.25, nudge_y = 0.25, 
                         check_overlap = T)
dev.off()

################################################################################################################################ 
################################################################################################################################ MDS plot using RLD
################################################################################################################################ 
################################################################################################################################

mds <- as.data.frame(colData(rld))  %>%  cbind(cmdscale(sampleDistMatrix))

png(paste(NAME, "part4.RLD.MDS.euclidean.png", sep="."))
ggplot(mds, aes(x = `1`, y = `2`, color = sample, shape = sample)) +
            geom_point(shape = c(1,2,3,4,5,6,7,8), size = 8) +
            coord_fixed() + 
            ggtitle("MDS with RLD data : euclidean") + 
            geom_text(    
                         label = samples$sample, 
                         nudge_x = 0.25, nudge_y = 0.25, 
                         check_overlap = T)
dev.off()

mdsPois <- as.data.frame(colData(dds)) %>% cbind(cmdscale(samplePoisDistMatrix))

png(paste(NAME, "part4.RLD.MDS.poisson.png", sep="."))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = sample, shape = sample)) +
            geom_point(shape = c(1,2,3,4,5,6,7,8), size = 8) +
            coord_fixed() + 
            ggtitle("MDS with RLD data : poisson") + 
            geom_text(    
                         label = samples$sample, 
                         nudge_x = 0.25, nudge_y = 0.25, 
                         check_overlap = T)
dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ DE ANALYSIS
################################################################################################################################
# it includes also the analysis with SUV
# it includes also the analysis with RUV-seq
# we may try the analysis with COMBAT-seq
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ Differentially Expressed Genes
################################################################################################################################ 
################################################################################################################################ 
# computing the DIFFERENTIAL EXPRESSION
# and extracting the CONTRASTS
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################

dds <- DESeq(dds)

resultsNames(dds)

# [1] "Intercept"              "group_lacSTZ_vs_lacCTR" "group_zfpCTR_vs_lacCTR"
# [4] "group_zfpSTZ_vs_lacCTR"

# using pre-existing size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

# res <- results(dds)

# head(res[order(res$pvalue),])
# tail(res[order(res$pvalue),])

# head(res[order(res$pvalue),])
# log2 fold change (MLE): group zfpSTZ vs lacCTR 
# Wald test p-value: group zfpSTZ vs lacCTR 
# DataFrame with 6 rows and 6 columns
#                     baseMean    log2FoldChange             lfcSE
#                    <numeric>         <numeric>         <numeric>
# Ddx3          1313.1303907189  9.02193651818142 0.836029219226851
# LOC100911498 1733.87019023218 -4.83960418567183 0.461054400769813
# Eif2s3y      2197.82733524275  14.5148361239418  1.49819308901979
# Fosl1        196.032857985014 -7.22399656906971  1.05013770635197
# Rgs1         169.474355049354  6.27315728776023   1.0047093895629
# Sectm1b      26.4652548079873  19.4095217311205  3.12976627475903
#                          stat               pvalue                 padj
#                     <numeric>            <numeric>            <numeric>
# Ddx3          10.7914129203819 3.77933841614821e-27 4.88252729982188e-23
# LOC100911498  -10.496818114286 8.93409710557449e-26 5.77098002534584e-22
# Eif2s3y       9.68822792624033 3.38348048440369e-22 1.45703947926704e-18
# tail(res[order(res$pvalue),])
# log2 fold change (MLE): group zfpSTZ vs lacCTR 
# Wald test p-value: group zfpSTZ vs lacCTR 
# DataFrame with 6 rows and 6 columns
#                     baseMean     log2FoldChange            lfcSE
#                    <numeric>          <numeric>        <numeric>
#Alyref       1257.94465509869   1.88067668847253 1.42882070064513
#Myo9b        97.6518631329639 -0.977596932574667 4.16821315911505
#LOC100911576 43.8462789493511   1.26356573503594 4.16956466830857
#Nek6         163.432885853561   2.38429486414968  2.9170146810241
#LOC108348108 913.166000542303  -2.08591632162596 2.76526410605146
#Chrna4       240.350229871719  0.162812123708916 1.04189897954528
#                           stat    pvalue      padj
#                      <numeric> <numeric> <numeric>
#Alyref         1.31624400991908        NA        NA
#Myo9b        -0.234536213781884        NA        NA
#LOC100911576  0.303045002429119        NA        NA
#Nek6          0.817374996313903        NA        NA
#LOC108348108 -0.754328064744765        NA        NA
#Chrna4        0.156264788530625        NA        NA

################################################################################################################################
################################################################################################################################
# The results table when printed will provide the information about
#     the comparison, e.g. "log2 fold change (MAP): condition treated vs
#     untreated", meaning that the estimates are of log2(treated /
#     untreated), as would be returned by
#     'contrast=c("condition","treated","untreated")'. Multiple results
#     can be returned for analyses beyond a simple two group comparison,
#     so 'results' takes arguments 'contrast' and 'name' to help the
#     user pick out the comparisons of interest for printing a results
#     table. The use of the 'contrast' argument is recommended for exact
#     specification of the levels which should be compared and their
#     order.

#     If 'results' is run without specifying 'contrast' or 'name', it
#     will return the comparison of the last level of the last variable
#     in the design formula over the first level of this variable. For
#     example, for a simple two-group comparison, this would return the
#     log2 fold changes of the second group over the first group (the
#     reference level). Please see examples below and in the vignette.

#     The argument 'contrast' can be used to generate results tables for
#     any comparison of interest, for example, the log2 fold change
#     between two levels of a factor, and its usage is described below.
#     It can also accomodate more complicated numeric comparisons. The
#     test statistic used for a contrast is:

#                        c' beta / sqrt( c' Sigma c )                    
#     
#     The argument 'name' can be used to generate results tables for
#     individual effects, which must be individual elements of
#     'resultsNames(object)'. These individual effects could represent
#     continuous covariates, effects for individual levels, or
#     individual interaction effects.


################################################################################################################################ plotCounts()
################################################################################################## to display SPECIFIC GENES : ZFP

png(paste(NAME, "part3.plotCounts.ZFP.png", sep="."))
plotCounts(dds, "ZFP", "group")
dev.off()

################################################################################################################################
################################################################################################## to search about plotMA ()
####################################################################### we will plot the RESULTS for the COMPARISON of INTEREST

# png(paste(NAME, "part3.MA.png", sep="."))
# plotMA(res, ylim=c(-5,5))
# dev.off()

# results(object, contrast, name, lfcThreshold = 0, ..)
# about the CONTRASTS : 

# log2 fold change (MLE): group zfpSTZ vs lacCTR 
# Wald test p-value: group zfpSTZ vs lacCTR 

################################################################################################################################
# resultsNames(dds)  ################################################################ showing the CONTRASTS that have been found 
# [1] "Intercept"              "group_lacSTZ_vs_lacCTR" "group_zfpCTR_vs_lacCTR"
# [4] "group_zfpSTZ_vs_lacCTR"
################################################################################################################################
#     results(object, contrast, name, lfcThreshold = 0,
#       altHypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
#       listValues = c(1, -1), cooksCutoff, independentFiltering = TRUE,
#       alpha = 0.1, filter, theta, pAdjustMethod = "BH", filterFun,
#       format = c("DataFrame", "GRanges", "GRangesList"), test,
#       addMLE = FALSE, tidy = FALSE, parallel = FALSE,
#       BPPARAM = bpparam(), minmu = 0.5)
     
#         resultsNames(object)
#         removeResults(object)
     
# object: a DESeqDataSet, on which one of the following functions has
#          already been called: 'DESeq', 'nbinomWaldTest', or
#          'nbinomLRT'

# contrast: this argument specifies what comparison to extract from the
#          'object' to build a results table. one of either:
#
#            • a character vector with exactly three elements: the name
#              of a factor in the design formula, the name of the
#              numerator level for the fold change, and the name of the
#              denominator level for the fold change (simplest case)
#
#            • a list of 2 character vectors: the names of the fold
#              changes for the numerator, and the names of the fold
#              changes for the denominator. these names should be
#              elements of 'resultsNames(object)'. if the list is length
#              1, a second element is added which is the empty character
#              vector, 'character()'. (more general case, can be to
#              combine interaction terms and main effects)
#
#            • a numeric contrast vector with one element for each
#              element in 'resultsNames(object)' (most general case)
#
#          If specified, the 'name' argument is ignored.
#
#    name: the name of the individual effect (coefficient) for building
#          a results table. Use this argument rather than 'contrast' for
#          continuous variables, individual effects or for individual
#          interaction terms. The value provided to 'name' must be an
#          element of 'resultsNames(object)'.
#    lfcThreshold: a non-negative value which specifies a log2 fold change
#          threshold. The default value is 0, corresponding to a test
#          that the log2 fold changes are equal to zero.

# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html

#[1] "Intercept"              "group_lacSTZ_vs_lacCTR" "group_zfpCTR_vs_lacCTR"
#[4] "group_zfpSTZ_vs_lacCTR"

# results(dds)

# plotCounts(dds, gene, intgroup = "condition", normalized = TRUE,
#       transform = TRUE, main, xlab = "group", returnData = FALSE,
#       replaced = FALSE, pc, ...)
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ IN ORDER TO EXTRACT the CONTRASTS of INTEREST :

# DESEq2 design : https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html

# in order to extract the CONTRASTS of INTEREST : FROM the GROUP, to extract the CONTRASTS of INTEREST 

# Calling the results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. 
# If there are more than 2 levels for this variable, results will extract the results table for a comparison of the last level over the first level. 

# unique(samples$group)
# [1] lacCTR lacSTZ zfpCTR zfpSTZ
# Levels: lacCTR lacSTZ zfpCTR zfpSTZ

# for example : 
# contrast: this argument specifies what comparison to extract from the
#          'object' to build a results table. one of either:
#
#            • a character vector with exactly three elements: the name
#              of a factor in the design formula, the name of the
#              numerator level for the fold change, and the name of the
#              denominator level for the fold change (simplest case)
# NUMERATOR
# DENOMINATOR

# lfcSE, the standard error estimate for the log2 fold change estimate. We can also express the uncertainty of a particular effect size estimate 
# as the result of a statistical test. The purpose of a test for differential expression is to test whether the data provides sufficient evidence 
# to conclude that this value is really different from zero. 

# results(dds, contrast=c("group", "zfpSTZ", "zfpCTR"))
# summary(res)

# out of 12927 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 20, 0.15%
# LFC < 0 (down)     : 6, 0.046%
# outliers [1]       : 8, 0.062%
# low counts [2]     : 0, 0%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

################################################################################################################################
################################################################################################################################ the CONTRASTS of INTEREST

lacCTR.vs.lacSTZ <- as.data.frame(results(dds, contrast=c("group", "lacSTZ", "lacCTR")))
head(lacCTR.vs.lacSTZ)
dim(lacCTR.vs.lacSTZ)
mcols(results(dds, contrast=c("group", "lacSTZ", "lacCTR")), use.names = TRUE)
lacCTR.vs.lacSTZ$GENE = rownames(lacCTR.vs.lacSTZ)

zfpCTR.vs.zfpSTZ <- as.data.frame(results(dds, contrast=c("group", "zfpSTZ", "zfpCTR")))
head(zfpCTR.vs.zfpSTZ)
dim(zfpCTR.vs.zfpSTZ)
mcols(results(dds, contrast=c("group", "zfpSTZ", "zfpCTR")), use.names = TRUE)
zfpCTR.vs.zfpSTZ$GENE = rownames(zfpCTR.vs.zfpSTZ)

lacCTR.vs.zfpCTR <- as.data.frame(results(dds, contrast=c("group", "zfpCTR", "lacCTR")))
head(lacCTR.vs.zfpCTR)
dim(lacCTR.vs.zfpCTR)
mcols(results(dds, contrast=c("group", "zfpCTR", "lacCTR")), use.names = TRUE)
lacCTR.vs.zfpCTR$GENE = rownames(lacCTR.vs.zfpCTR)

lacSTZ.vs.zfpSTZ <- as.data.frame(results(dds, contrast=c("group", "zfpSTZ", "lacSTZ")))
head(lacSTZ.vs.zfpSTZ)
dim(lacSTZ.vs.zfpSTZ)
mcols(results(dds, contrast=c("group", "zfpSTZ", "lacSTZ")), use.names = TRUE)
lacSTZ.vs.zfpSTZ$GENE = rownames(lacSTZ.vs.zfpSTZ)

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ to print the RESULTS in FILES 

write.table(as.data.frame(lacCTR.vs.lacSTZ),
            file=paste(NAME, "differential.expression.lacCTR.vs.lacSTZ", "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(zfpCTR.vs.zfpSTZ),
            file=paste(NAME, "differential.expression.zfpCTR.vs.zfpSTZ", "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacCTR.vs.zfpCTR),
            file=paste(NAME, "differential.expression.lacCTR.vs.zfpCTR", "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacSTZ.vs.zfpSTZ),
            file=paste(NAME, "differential.expression.lacSTZ.vs.zfpSTZ", "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ COMPUTING the NUMBER of the REGULATED GENES 

colnames(lacCTR.vs.lacSTZ)
colnames(zfpCTR.vs.zfpSTZ)
colnames(lacCTR.vs.zfpCTR)
colnames(lacSTZ.vs.zfpSTZ)

dim(subset(lacCTR.vs.lacSTZ, (lacCTR.vs.lacSTZ.pvalue < 0.05) & (abs(lacCTR.vs.lacSTZ.log2FoldChange) > 2)))
dim(subset(lacCTR.vs.lacSTZ, (lacCTR.vs.lacSTZ.padj < 0.1)  & (abs(lacCTR.vs.lacSTZ.log2FoldChange) > 2)))

dim(subset(zfpCTR.vs.zfpSTZ, (zfpCTR.vs.zfpSTZ.pvalue < 0.05)  & (abs(zfpCTR.vs.zfpSTZ.log2FoldChange) > 2)))
dim(subset(zfpCTR.vs.zfpSTZ, (zfpCTR.vs.zfpSTZ.padj < 0.1)  & (abs(zfpCTR.vs.zfpSTZ.log2FoldChange) > 2)))
 
dim(subset(lacCTR.vs.zfpCTR, (lacCTR.vs.zfpCTR.pvalue < 0.05)  & (abs(lacCTR.vs.zfpCTR.log2FoldChange) > 2)))
dim(subset(lacCTR.vs.zfpCTR, (lacCTR.vs.zfpCTR.padj < 0.1)  & (abs(lacCTR.vs.zfpCTR.log2FoldChange) > 2)))

dim(subset(lacSTZ.vs.zfpSTZ,  (lacSTZ.vs.zfpSTZ.pvalue < 0.05)  & (abs(lacSTZ.vs.zfpSTZ.log2FoldChange) > 2)))
dim(subset(lacSTZ.vs.zfpSTZ,  (lacSTZ.vs.zfpSTZ.padj < 0.1)  & (abs(lacSTZ.vs.zfpSTZ.log2FoldChange) > 2)))

dim(subset(lacCTR.vs.lacSTZ, (pvalue < 0.05) & (abs(log2FoldChange) > 1.2)))
dim(subset(lacCTR.vs.lacSTZ, (padj < 0.1)  & (abs(log2FoldChange) > 1.2)))

dim(subset(zfpCTR.vs.zfpSTZ, (pvalue < 0.05)  & (abs(log2FoldChange) > 1.2)))
dim(subset(zfpCTR.vs.zfpSTZ, (padj < 0.1)  & (abs(log2FoldChange) > 1.2)))
 
dim(subset(lacCTR.vs.zfpCTR, (pvalue < 0.05)  & (abs(log2FoldChange) > 1.2)))
dim(subset(lacCTR.vs.zfpCTR, (padj < 0.1)  & (abs(log2FoldChange) > 1.2)))

dim(subset(lacSTZ.vs.zfpSTZ,  (pvalue < 0.05)  & (abs(log2FoldChange) > 1.2)))
dim(subset(lacSTZ.vs.zfpSTZ,  (padj < 0.1)  & (abs(log2FoldChange) > 1.2)))

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ the MA plots
################################################################################################################################ 

resultsNames(dds)
resultsNames(dds)

lacCTR.vs.lacSTZ.lfc <- as.data.frame(lfcShrink(dds, contrast=c("group", "lacSTZ", "lacCTR")))
lacCTR.vs.lacSTZ.lfc$GENE = rownames(lacCTR.vs.lacSTZ.lfc)

zfpCTR.vs.zfpSTZ.lfc <- as.data.frame(lfcShrink(dds, contrast=c("group", "zfpSTZ", "zfpCTR")))
zfpCTR.vs.zfpSTZ.lfc$GENE = rownames(zfpCTR.vs.zfpSTZ.lfc)

lacCTR.vs.zfpCTR.lfc <- as.data.frame(lfcShrink(dds, contrast=c("group", "zfpCTR", "lacCTR")))
lacCTR.vs.zfpCTR.lfc$GENE = rownames(lacCTR.vs.zfpCTR.lfc)

lacSTZ.vs.zfpSTZ.lfc <- as.data.frame(lfcShrink(dds, contrast=c("group", "zfpSTZ", "lacSTZ")))
lacSTZ.vs.zfpSTZ.lfc$GENE = rownames(lacSTZ.vs.zfpSTZ.lfc)

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 

dim(subset(lacCTR.vs.lacSTZ.lfc, pvalue < 0.05))
dim(subset(lacCTR.vs.lacSTZ.lfc, padj < 0.1))

dim(subset(zfpCTR.vs.zfpSTZ.lfc, pvalue < 0.05))
dim(subset(zfpCTR.vs.zfpSTZ.lfc, padj < 0.1))
 
dim(subset(lacCTR.vs.zfpCTR.lfc, pvalue < 0.05))
dim(subset(lacCTR.vs.zfpCTR.lfc, padj < 0.1))

dim(subset(lacSTZ.vs.zfpSTZ.lfc,  pvalue < 0.05))
dim(subset(lacSTZ.vs.zfpSTZ.lfc, padj < 0.1))

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################

write.table(as.data.frame(lacCTR.vs.lacSTZ.lfc),
            file=paste(NAME, "differential.expression.lacCTR.vs.lacSTZ", "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(zfpCTR.vs.zfpSTZ.lfc),
            file=paste(NAME, "differential.expression.zfpCTR.vs.zfpSTZ", "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacCTR.vs.zfpCTR.lfc),
            file=paste(NAME, "differential.expression.lacCTR.vs.zfpCTR", "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacSTZ.vs.zfpSTZ.lfc),
            file=paste(NAME, "differential.expression.lacSTZ.vs.zfpSTZ", "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ the PIECE of CODE is valid

#PATH = getwd()
#library("ReportingTools")
#res = results(dds, contrast=c("group", "lacSTZ", "lacCTR"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here
#rep <- HTMLReport(shortName = paste(NAME, "lacCTR.vs.lacSTZ", sep="."),  
#                  title = ,
#                  basePath = tmp, 
#                  reportDirectory = "report")
#publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
#finish(rep)

# it launches the report in a web browser:
# browseURL(file.path(tmp,"report",".html"))

################################################################################################################################
################################################################################################################################

#PATH = getwd()
#library("ReportingTools")
#res = results(dds, contrast=c("group", "zfpSTZ", "zfpCTR"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here
#rep <- HTMLReport(shortName = paste(NAME, "zfpCTR.vs.zfpSTZ", sep="."),  
#                  title = ,
#                  basePath = tmp, 
#                  reportDirectory = "report")
#publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
#finish(rep)

# it launches the report in a web browser:
# browseURL(file.path(tmp,"report",".html"))

################################################################################################################################
################################################################################################################################

#PATH = getwd()
#library("ReportingTools")
#res = results(dds, contrast=c("group", "zfpCTR", "lacCTR"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here
#rep <- HTMLReport(shortName = paste(NAME, "lacCTR.vs.zfpCTR", sep="."),  
#                  title = ,
#                  basePath = tmp, 
#                  reportDirectory = "report")
#publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
#finish(rep)

# it launches the report in a web browser:
# browseURL(file.path(tmp,"report",".html"))

################################################################################################################################
################################################################################################################################

#PATH = getwd()
#library("ReportingTools")
#res = results(dds, contrast=c("group", "zfpSTZ", "lacSTZ"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here
#rep <- HTMLReport(shortName = paste(NAME, "lacSTZ.vs.zfpSTZ", sep="."),  
#                  title = ,
#                  basePath = tmp, 
#                  reportDirectory = "report")
#publish(res, rep, dds, n=500, make.plots=TRUE, factor=dds$group)
#finish(rep)

# it launches the report in a web browser:
# browseURL(file.path(tmp,"report",".html"))

################################################################################################################################
################################################################################################################################

### the RESULTS are in the FOLDER :
###"/tmp/report/KALLISTO.RefSeq.INTEGRATED.TPM.txt.no.dup.genes.txt.lacSTZ.vs.zfpSTZ.html"
 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ writing the REPORTS

#PATH = getwd()
#library("Glimma")
#res = results(dds, contrast=c("group", "lacSTZ", "lacCTR"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here

#status <- as.numeric(res$padj < 0.3)
#anno <- data.frame(GeneID=rownames(res), symbol=rownames(res))
#glMDPlot(res, 
#         status=status, 
#         counts=counts(dds,normalized=TRUE),
#         groups=dds$group, 
#         transform=FALSE,
#         samples=colnames(dds), 
#         anno=anno,
#         path=tmp, 
#         folder="glimma", 
#         launch=FALSE)

# to launch the report in a web browser:
# browseURL(file.path(tmp,"glimma","MD-Plot.html"))

#PATH = getwd()
#library("Glimma")
#res = results(dds, contrast=c("group", "zfpSTZ", "zfpCTR"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here

#status <- as.numeric(res$padj < 0.3)
#anno <- data.frame(GeneID=rownames(res), symbol=rownames(res))
#glMDPlot(res, 
#         status=status, 
#         counts=counts(dds,normalized=TRUE),
#         groups=dds$group, 
#         transform=FALSE,
#         samples=colnames(dds), 
#         anno=anno,
#         path=tmp, 
#         folder="glimma", 
#         launch=FALSE)

# to launch the report in a web browser:
# browseURL(file.path(tmp,"glimma","MD-Plot.html"))

#PATH = getwd()
#library("Glimma")
#res = results(dds, contrast=c("group", "zfpCTR", "lacCTR"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here

#status <- as.numeric(res$padj < 0.3)
#anno <- data.frame(GeneID=rownames(res), symbol=rownames(res))
#glMDPlot(res, 
#         status=status, 
#         counts=counts(dds,normalized=TRUE),
#         groups=dds$group, 
#         transform=FALSE,
#         samples=colnames(dds), 
#         anno=anno,
#         path=tmp, 
#         folder="glimma", 
#         launch=FALSE)

# to launch the report in a web browser:
# browseURL(file.path(tmp,"glimma","MD-Plot.html"))

#PATH = getwd()
#library("Glimma")
#res = results(dds, contrast=c("group", "zfpSTZ", "lacSTZ"))
#tmp <- tempdir(PATH) # you would instead use a meaningful path here

#status <- as.numeric(res$padj < 0.3)
#anno <- data.frame(GeneID=rownames(res), symbol=rownames(res))
#glMDPlot(res, 
#         status=status, 
#         counts=counts(dds,normalized=TRUE),
#         groups=dds$group, 
#         transform=FALSE,
#         samples=colnames(dds), 
#         anno=anno,
#         path=tmp, 
#         folder="glimma", 
#         launch=FALSE)

# to launch the report in a web browser:
# browseURL(file.path(tmp,"glimma","MD-Plot.html"))

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################
################################################################################################################################ INTEGRATING these DATAFRAMES
################################################################################################################################
################################################################################################################################ INTEGRATING these FILES 

# using the COUNTS
# counts = read.delim("KALLISTO.RefSeq.INTEGRATED.TPM.txt.no.dup.genes.txt", 
#                      sep="\t", header=T, stringsAsFactors=F)

# lacCTR.vs.lacSTZ <- as.data.frame(results(dds, contrast=c("group", "lacSTZ", "lacCTR")))
# zfpCTR.vs.zfpSTZ <- as.data.frame(results(dds, contrast=c("group", "zfpSTZ", "zfpCTR")))
# lacCTR.vs.zfpCTR <- as.data.frame(results(dds, contrast=c("group", "zfpCTR", "lacCTR")))
# lacSTZ.vs.zfpSTZ <- as.data.frame(results(dds, contrast=c("group", "zfpSTZ", "lacSTZ")))

# using these DATAFRAMES

# using COUNTS
# lacCTR.vs.lacSTZ 
# zfpCTR.vs.zfpSTZ 
# lacCTR.vs.zfpCTR 
# lacSTZ.vs.zfpSTZ 

dim(counts)
dim(lacCTR.vs.lacSTZ) 
dim(zfpCTR.vs.zfpSTZ) 
dim(lacCTR.vs.zfpCTR) 
dim(lacSTZ.vs.zfpSTZ) 

png(paste(NAME, "part3.MA", "lacCTR.vs.lacSTZ", "png", sep="."))
plot(log2(lacCTR.vs.lacSTZ$baseMean), lacCTR.vs.lacSTZ$log2FoldChange)
lines(lowess(log2(lacCTR.vs.lacSTZ$baseMean), lacCTR.vs.lacSTZ$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

png(paste(NAME, "part3.MA", "zfpCTR.vs.zfpSTZ", "png", sep="."))
plot(log2(zfpCTR.vs.zfpSTZ$baseMean), zfpCTR.vs.zfpSTZ$log2FoldChange)
lines(lowess(log2(zfpCTR.vs.zfpSTZ$baseMean), zfpCTR.vs.zfpSTZ$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

png(paste(NAME, "part3.MA", "lacCTR.vs.zfpCTR", "png", sep="."))
plot(log2(lacCTR.vs.zfpCTR$baseMean), lacCTR.vs.zfpCTR$log2FoldChange)
lines(lowess(log2(lacCTR.vs.zfpCTR$baseMean), lacCTR.vs.zfpCTR$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

png(paste(NAME, "part3.MA", "lacSTZ.vs.zfpSTZ", "png", sep="."))
plot(log2(lacSTZ.vs.zfpSTZ$baseMean), lacSTZ.vs.zfpSTZ$log2FoldChange)
lines(lowess(log2(lacSTZ.vs.zfpSTZ$baseMean), lacSTZ.vs.zfpSTZ$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

colnames(counts)
colnames(lacCTR.vs.lacSTZ) 
colnames(zfpCTR.vs.zfpSTZ) 
colnames(lacCTR.vs.zfpCTR) 
colnames(lacSTZ.vs.zfpSTZ) 

# colnames(counts)
# [1] "target_id.S1"  "length.S1"     "eff_length.S1" "est_counts.S1"
# [5] "tpm.S1"        "SYMBOL.x"      "length.S2"     "eff_length.S2"
# [9] "est_counts.S2" "tpm.S2"        "SYMBOL.y"      "length.S3"    
#[13] "eff_length.S3" "est_counts.S3" "tpm.S3"        "SYMBOL.x.1"   
#[17] "length.S4"     "eff_length.S4" "est_counts.S4" "tpm.S4"       
#[21] "SYMBOL.y.1"    "length.S5"     "eff_length.S5" "est_counts.S5"
#[25] "tpm.S5"        "SYMBOL.x.2"    "length.S6"     "eff_length.S6"
#[29] "est_counts.S6" "tpm.S6"        "SYMBOL.y.2"    "length.S7"    
#[33] "eff_length.S7" "est_counts.S7" "tpm.S7"        "SYMBOL.x.3"   
#[37] "length.S8"     "eff_length.S8" "est_counts.S8" "tpm.S8"       
#[41] "SYMBOL.y.3"    "length.S9"     "eff_length.S9" "est_counts.S9"
#[45] "tpm.S9"        "SYMBOL"       
# colnames(lacCTR.vs.lacSTZ) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE"          
# colnames(zfpCTR.vs.zfpSTZ) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE"          
# colnames(lacCTR.vs.zfpCTR) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE"          
# colnames(lacSTZ.vs.zfpSTZ) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE" 

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
# using COUNTS
# lacCTR.vs.lacSTZ 
# zfpCTR.vs.zfpSTZ 
# lacCTR.vs.zfpCTR 
# lacSTZ.vs.zfpSTZ 
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

counts.SMALL.TPM = counts.TPM ####################################################### to use this DATAFRAME for DATA INTEGRATION 

dim(counts.SMALL.TPM)
head(counts.SMALL.TPM)
tail(counts.SMALL.TPM)

dim(lacCTR.vs.lacSTZ) 
dim(zfpCTR.vs.zfpSTZ) 
dim(lacCTR.vs.zfpCTR) 
dim(lacSTZ.vs.zfpSTZ) 

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

# baseMean log2FoldChange     lfcSE       stat     pvalue      padj

# dim(lacCTR.vs.lacSTZ) 
# dim(zfpCTR.vs.zfpSTZ) 
# dim(lacCTR.vs.zfpCTR) 
# dim(lacSTZ.vs.zfpSTZ) 

head(lacCTR.vs.lacSTZ) 

colnames(lacCTR.vs.lacSTZ)[1] = "lacCTR.vs.lacSTZ.baseMean"
colnames(lacCTR.vs.lacSTZ)[2] = "lacCTR.vs.lacSTZ.log2FoldChange"
colnames(lacCTR.vs.lacSTZ)[3] = "lacCTR.vs.lacSTZ.lfcSE"
colnames(lacCTR.vs.lacSTZ)[4] = "lacCTR.vs.lacSTZ.stat"
colnames(lacCTR.vs.lacSTZ)[5] = "lacCTR.vs.lacSTZ.pvalue"
colnames(lacCTR.vs.lacSTZ)[6] = "lacCTR.vs.lacSTZ.padj"

head(zfpCTR.vs.zfpSTZ) 

colnames(zfpCTR.vs.zfpSTZ)[1] = "zfpCTR.vs.zfpSTZ.baseMean"
colnames(zfpCTR.vs.zfpSTZ)[2] = "zfpCTR.vs.zfpSTZ.log2FoldChange"
colnames(zfpCTR.vs.zfpSTZ)[3] = "zfpCTR.vs.zfpSTZ.lfcSE"
colnames(zfpCTR.vs.zfpSTZ)[4] = "zfpCTR.vs.zfpSTZ.stat"
colnames(zfpCTR.vs.zfpSTZ)[5] = "zfpCTR.vs.zfpSTZ.pvalue"
colnames(zfpCTR.vs.zfpSTZ)[6] = "zfpCTR.vs.zfpSTZ.padj"

head(lacCTR.vs.zfpCTR) 

colnames(lacCTR.vs.zfpCTR)[1] = "lacCTR.vs.zfpCTR.baseMean"
colnames(lacCTR.vs.zfpCTR)[2] = "lacCTR.vs.zfpCTR.log2FoldChange"
colnames(lacCTR.vs.zfpCTR)[3] = "lacCTR.vs.zfpCTR.lfcSE"
colnames(lacCTR.vs.zfpCTR)[4] = "lacCTR.vs.zfpCTR.stat"
colnames(lacCTR.vs.zfpCTR)[5] = "lacCTR.vs.zfpCTR.pvalue"
colnames(lacCTR.vs.zfpCTR)[6] = "lacCTR.vs.zfpCTR.padj"

head(lacSTZ.vs.zfpSTZ) 

colnames(lacSTZ.vs.zfpSTZ)[1] = "lacSTZ.vs.zfpSTZ.baseMean"
colnames(lacSTZ.vs.zfpSTZ)[2] = "lacSTZ.vs.zfpSTZ.log2FoldChange"
colnames(lacSTZ.vs.zfpSTZ)[3] = "lacSTZ.vs.zfpSTZ.lfcSE"
colnames(lacSTZ.vs.zfpSTZ)[4] = "lacSTZ.vs.zfpSTZ.stat"
colnames(lacSTZ.vs.zfpSTZ)[5] = "lacSTZ.vs.zfpSTZ.pvalue"
colnames(lacSTZ.vs.zfpSTZ)[6] = "lacSTZ.vs.zfpSTZ.padj"

###############################################################
###############################################################
###############################################################

head(lacCTR.vs.lacSTZ) 
head(zfpCTR.vs.zfpSTZ) 
head(lacCTR.vs.zfpCTR) 
head(lacSTZ.vs.zfpSTZ) 

###############################################################
############################################################### we have to do these again
###############################################################

counts.TPM <- subset(counts, select=c("SYMBOL", 
                "tpm.S1",
                "tpm.S2",
                "tpm.S3", 
                "tpm.S4", 
                "tpm.S5",
                "tpm.S6", 
                "tpm.S7", 
                "tpm.S8"))

head(counts.TPM)
dim(counts.TPM)  

rownames(counts.TPM) = counts.TPM$SYMBOL

################################################################################################################################ 
################################################################################################################################
######################################################## to use the following for the DATA INTEGRATION at the end of the script:

counts.SMALL.TPM = counts.TPM ####################################################### to use this DATAFRAME for DATA INTEGRATION 

dim(counts.SMALL.TPM)

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

df1 = merge(counts.SMALL.TPM, 
            lacCTR.vs.lacSTZ, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

df2 = merge(df1, 
            zfpCTR.vs.zfpSTZ, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

df3 = merge(df2,
            lacCTR.vs.zfpCTR, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

df4 = merge(df3,
            lacSTZ.vs.zfpSTZ, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

write.table(df4,
            file=paste(NAME, "the.INTEGRATED.TABLES", "RESULTS.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

# baseMean log2FoldChange     lfcSE       stat     pvalue      padj

dim(lacCTR.vs.lacSTZ.lfc) 
dim(zfpCTR.vs.zfpSTZ.lfc) 
dim(lacCTR.vs.zfpCTR.lfc) 
dim(lacSTZ.vs.zfpSTZ.lfc) 

################################################################################################################################
################################################################################################################################

head(lacCTR.vs.lacSTZ.lfc) 

colnames(lacCTR.vs.lacSTZ.lfc)[1] = "lacCTR.vs.lacSTZ.baseMean.lfc"
colnames(lacCTR.vs.lacSTZ.lfc)[2] = "lacCTR.vs.lacSTZ.log2FoldChange.lfc"
colnames(lacCTR.vs.lacSTZ.lfc)[3] = "lacCTR.vs.lacSTZ.lfcSE.lfc"
colnames(lacCTR.vs.lacSTZ.lfc)[4] = "lacCTR.vs.lacSTZ.stat.lfc"
colnames(lacCTR.vs.lacSTZ.lfc)[5] = "lacCTR.vs.lacSTZ.pvalue.lfc"
colnames(lacCTR.vs.lacSTZ.lfc)[6] = "lacCTR.vs.lacSTZ.padj.lfc"

head(zfpCTR.vs.zfpSTZ.lfc) 

colnames(zfpCTR.vs.zfpSTZ.lfc)[1] = "zfpCTR.vs.zfpSTZ.baseMean.lfc"
colnames(zfpCTR.vs.zfpSTZ.lfc)[2] = "zfpCTR.vs.zfpSTZ.log2FoldChange.lfc"
colnames(zfpCTR.vs.zfpSTZ.lfc)[3] = "zfpCTR.vs.zfpSTZ.lfcSE.lfc"
colnames(zfpCTR.vs.zfpSTZ.lfc)[4] = "zfpCTR.vs.zfpSTZ.stat.lfc"
colnames(zfpCTR.vs.zfpSTZ.lfc)[5] = "zfpCTR.vs.zfpSTZ.pvalue.lfc"
colnames(zfpCTR.vs.zfpSTZ.lfc)[6] = "zfpCTR.vs.zfpSTZ.padj.lfc"

head(lacCTR.vs.zfpCTR.lfc) 

colnames(lacCTR.vs.zfpCTR.lfc)[1] = "lacCTR.vs.zfpCTR.baseMean.lfc"
colnames(lacCTR.vs.zfpCTR.lfc)[2] = "lacCTR.vs.zfpCTR.log2FoldChange.lfc"
colnames(lacCTR.vs.zfpCTR.lfc)[3] = "lacCTR.vs.zfpCTR.lfcSE.lfc"
colnames(lacCTR.vs.zfpCTR.lfc)[4] = "lacCTR.vs.zfpCTR.stat.lfc"
colnames(lacCTR.vs.zfpCTR.lfc)[5] = "lacCTR.vs.zfpCTR.pvalue.lfc"
colnames(lacCTR.vs.zfpCTR.lfc)[6] = "lacCTR.vs.zfpCTR.padj.lfc"

head(lacSTZ.vs.zfpSTZ.lfc) 

colnames(lacSTZ.vs.zfpSTZ.lfc)[1] = "lacSTZ.vs.zfpSTZ.baseMean.lfc"
colnames(lacSTZ.vs.zfpSTZ.lfc)[2] = "lacSTZ.vs.zfpSTZ.log2FoldChange.lfc"
colnames(lacSTZ.vs.zfpSTZ.lfc)[3] = "lacSTZ.vs.zfpSTZ.lfcSE.lfc"
colnames(lacSTZ.vs.zfpSTZ.lfc)[4] = "lacSTZ.vs.zfpSTZ.stat.lfc"
colnames(lacSTZ.vs.zfpSTZ.lfc)[5] = "lacSTZ.vs.zfpSTZ.pvalue.lfc"
colnames(lacSTZ.vs.zfpSTZ.lfc)[6] = "lacSTZ.vs.zfpSTZ.padj.lfc"

###############################################################
###############################################################

head(lacCTR.vs.lacSTZ.lfc) 
head(zfpCTR.vs.zfpSTZ.lfc) 
head(lacCTR.vs.zfpCTR.lfc) 
head(lacSTZ.vs.zfpSTZ.lfc) 

###############################################################
###############################################################

df1.lfc = merge(counts.SMALL.TPM, 
            lacCTR.vs.lacSTZ.lfc, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

df2.lfc = merge(df1.lfc, 
            zfpCTR.vs.zfpSTZ.lfc, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

df3.lfc = merge(df2.lfc,
            lacCTR.vs.zfpCTR.lfc, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

df4.lfc = merge(df3.lfc,
            lacSTZ.vs.zfpSTZ.lfc, 
            by.x = "SYMBOL", 
            by.y = "GENE", 
            all.x = TRUE)

write.table(df4.lfc,
            file=paste(NAME, "the.INTEGRATED.TABLES", "RESULTS.LFC.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

###############################################################
###############################################################
###############################################################

head(df1.lfc) 
head(df2.lfc) 
head(df3.lfc) 
head(df4.lfc) 

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
library(sva)
library(RUVSeq)
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ using SVAseq

dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ group, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))

svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

########################################################################################################################
########################################################################################################################

png(paste(NAME, "part5.SVA.png", sep="."))
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$group, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
dev.off()

# Finally, in order to use SVA to remove any effect on the counts from our surrogate variables, we simply add these two surrogate variables 
# as columns to the DESeqDataSet and then add them to the design:

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + group
  
ddssva$SV1
ddssva$SV2

# length(ddssva$SV1)
# length(ddssva$SV2)

########################################################################################################################
########################################################################################################################

ddssva 
colData(ddssva) 
assay(ddssva) 
design(ddssva) 

######################################## in order to make the transition to the next STEP :

ddssva <- DESeq(ddssva)

resultsNames(ddssva)

# vsd2 <- vst(ddssva, blind=TRUE)
# plotPCA(vsd2, "group")

# rld2 <- rlog(ddssva, blind = TRUE)
# plotPCA(rld2, "group")

########################################################################################################################

source("script_analysis_SVA.R")

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ using RUV-seq
# library("RUVSeq")

# we pull out a set of empirical control genes by looking at the genes that do not have a small p-value.

# set <- newSeqExpressionSet(counts(dds))
# idx  <- rowSums(counts(set) > 5) >= 2
# set  <- set[idx, ]
# set <- betweenLaneNormalization(set, which="upper")
# not.sig <- rownames(res)[which(res$pvalue > .1)]
# empirical <- rownames(set)[ rownames(set) %in% not.sig ]
# set <- RUVg(set, empirical, k=2)
# pData(set)

# pData(set)
#           W_1         W_2
#S1 -0.16744408 -0.29075436
#S2  0.61103370  0.14817442
#S3  0.11717069  0.01607802
#S4 -0.02289371  0.10000707
#S5 -0.34762428  0.09112168
#S6 -0.46130851  0.54105598
#S7 -0.24037062 -0.26050671
#S8  0.07945553 -0.64943424
#S9  0.43198127  0.30425814

#png(paste(NAME, "part5.RUV.png", sep="."))
#par(mfrow = c(2, 1), mar = c(3,5,3,1))
#for (i in 1:2) {
#  stripchart(pData(set)[, i] ~ dds$group, vertical = TRUE, main = paste0("W", i))
#  abline(h = 0)
#}
#dev.off()

# As before, if we wanted to control for these factors, we simply add them to the DESeqDataSet and to the design:

#ddsruv <- dds
#ddsruv$W1 <- set$W_1
#ddsruv$W2 <- set$W_2
#design(ddsruv) <- ~ W1 + W2 + group

########################################################################################################################
########################################################################################################################
########################################################################################################################

#ddsruv 
#colData(ddsruv) 
#assay(ddsruv) 
#design(ddsruv) 

########################################################################################################################
########################################################################################################################
########################################################################################################################

#ddsruv <- DESeq(ddsruv)

# vsd3 <- vst(ddsruv, blind=TRUE)
# plotPCA(vsd3, "group")

# rld3 <- rlog(ddsruv, blind = TRUE)
# plotPCA(rld3, "group")

########################################################################################################################

# source("script_analysis_RUV.R")

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
