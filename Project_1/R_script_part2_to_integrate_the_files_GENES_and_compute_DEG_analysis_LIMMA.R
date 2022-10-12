############################################################################################
############################################################################################ 
### it was a previous version (v7) of the script

library("ggplot2")
library("reshape2")
library("data.table")
library("tidyr")
library("gplots")
library("pheatmap")
library("RColorBrewer")

library("limma")
library("edgeR")
library("Glimma")
library("DESeq2")
library(data.table)

############################################################################################
############################################################################################
######################################## reading the files with the GENE EXPRESSION COUNTS :

genes <- read.delim("the_GENES.58381_genes.gencode.v28.basic.annotation.28aug2018.txt.INTEGRATED.file.ALL.samples.ALL.fields.txt",
                     sep="\t", header=T, stringsAsFactors=F)

head(genes) 
dim(genes)

############################################################################################
############################################################################################
############################################ transforming the DATA FRAME into a DATA TABLE :

genes.dt <- as.data.table(genes)

head(genes.dt) 
dim(genes.dt)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
##### in the next section, we are going to select the COLUMNS with COUNTS for DEG analysis :

# colnames(genes) 
# [1] "CHR"               "START"             "END"              
# [4] "STRAND"            "GENE_ID"           "GENE_NAME"        
# [7] "GENE_TYPE"         "DMSO1_lane1.count" "DMSO1_lane1.TPM"  
#[10] "DMSO1_lane1.FPKM"  "DMSO1_lane2.count" "DMSO1_lane2.TPM"  
#[13] "DMSO1_lane2.FPKM"  "DMSO2_lane1.count" "DMSO2_lane1.TPM"  
#[16] "DMSO2_lane1.FPKM"  "DMSO2_lane2.count" "DMSO2_lane2.TPM"  
#[19] "DMSO2_lane2.FPKM"  "DMSO3_lane1.count" "DMSO3_lane1.TPM"  
#[22] "DMSO3_lane1.FPKM"  "DMSO3_lane2.count" "DMSO3_lane2.TPM"  
#[25] "DMSO3_lane2.FPKM"  "Aph1.count"        "Aph1.TPM"         
#[28] "Aph1.FPKM"         "Aph2.count"        "Aph2.TPM"         
#[31] "Aph2.FPKM"         "Aph3.count"        "Aph3.TPM"         
#[34] "Aph3.FPKM"         "Aph_KH7_1.count"   "Aph_KH7_1.TPM"    
#[37] "Aph_KH7_1.FPKM"    "Aph_KH7_2.count"   "Aph_KH7_2.TPM"    
#[40] "Aph_KH7_2.FPKM"    "Aph_KH7_3.count"   "Aph_KH7_3.TPM"    
#[43] "Aph_KH7_3.FPKM"    "KH7_1.count"       "KH7_1.TPM"        
#[46] "KH7_1.FPKM"        "KH7_2.count"       "KH7_2.TPM"        
#[49] "KH7_2.FPKM"        "KH7_3.count"       "KH7_3.TPM"        
#[52] "KH7_3.FPKM"        "Noc_1.count"       "Noc_1.TPM"        
#[55] "Noc_1.FPKM"        "Noc_2.count"       "Noc_2.TPM"        
#[58] "Noc_2.FPKM"        "Noc_3.count"       "Noc_3.TPM"        
#[61] "Noc_3.FPKM"       

############################################################################################
############################################################################################

genes$ID <- rownames(genes) 
genes$GENE_NAME_ID <- paste(genes$GENE_NAME, 
                            genes$ID, sep=":")

head(genes)
dim(genes)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
####################################################### making a DATAFRAME of GENES COUNTS :

genes.counts <- subset(genes, select=c("GENE_NAME_ID", 
                       "DMSO1_lane1.count", "DMSO1_lane2.count",
                       "DMSO2_lane1.count", "DMSO2_lane2.count",
                       "DMSO3_lane1.count", "DMSO3_lane2.count", 
                       "Aph1.count", "Aph2.count", "Aph3.count", 
                       "Aph_KH7_1.count","Aph_KH7_2.count","Aph_KH7_3.count",
                       "KH7_1.count", "KH7_2.count", "KH7_3.count",
                       "Noc_1.count", "Noc_2.count", "Noc_3.count" ))

rownames(genes.counts) <- genes.counts$GENE_NAME_ID 
genes.counts <- genes.counts[,-1]

head(genes.counts)
dim(genes.counts)

############################################################################################
############################################################################################
############################################################################################
########################################################## making a DATAFRAME based on TPM :

genes.tpm <- subset(genes, select=c("GENE_NAME_ID", 
                       "DMSO1_lane1.TPM", "DMSO1_lane2.TPM",
                       "DMSO2_lane1.TPM", "DMSO2_lane2.TPM",
                       "DMSO3_lane1.TPM", "DMSO3_lane2.TPM", 
                       "Aph1.TPM", "Aph2.TPM", "Aph3.TPM", 
                       "Aph_KH7_1.TPM","Aph_KH7_2.TPM","Aph_KH7_3.TPM",
                       "KH7_1.TPM", "KH7_2.TPM", "KH7_3.TPM",
                       "Noc_1.TPM", "Noc_2.TPM", "Noc_3.TPM" ))

rownames(genes.tpm) <- genes.tpm$GENE_NAME_ID 
genes.tpm <- genes.tpm[,-1]

head(genes.tpm)
dim(genes.tpm)

############################################################################################
############################################################################################
############################################################################################
######################################################### making a DATAFRAME based on FPKM :

genes.fpkm <- subset(genes, select=c("GENE_NAME_ID", 
                       "DMSO1_lane1.FPKM", "DMSO1_lane2.FPKM",
                       "DMSO2_lane1.FPKM", "DMSO2_lane2.FPKM",
                       "DMSO3_lane1.FPKM", "DMSO3_lane2.FPKM", 
                       "Aph1.FPKM", "Aph2.FPKM", "Aph3.FPKM", 
                       "Aph_KH7_1.FPKM","Aph_KH7_2.FPKM","Aph_KH7_3.FPKM",
                       "KH7_1.FPKM", "KH7_2.FPKM", "KH7_3.FPKM",
                       "Noc_1.FPKM", "Noc_2.FPKM", "Noc_3.FPKM" ))

rownames(genes.fpkm) <- genes.fpkm$GENE_NAME_ID 
genes.fpkm <- genes.fpkm[,-1]

head(genes.fpkm)
dim(genes.fpkm)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
################# continuing to work with the DATAFRAME containing the COUNTS : genes.counts
################# in order to assess the DIFFERENTIAL EXPRESSION

head(genes.counts)
dim(genes.counts)

############################################################################################
############################################################################################
############################################################################################
######################################################## and SUBSETTING by SPECIFIC SAMPLES :

genes.counts.Aph <- subset(genes.counts, select=c(
                           "DMSO1_lane1.count", 
                           "DMSO2_lane1.count", 
                           "DMSO3_lane1.count",  
                           "Aph1.count", "Aph2.count", "Aph3.count" ))

dim(genes.counts.Aph)
head(genes.counts.Aph)

###############################################################################
###############################################################################

genes.counts.Aph_KH7 <- subset(genes.counts, select=c( 
                               "DMSO1_lane1.count", 
                               "DMSO2_lane1.count", 
                               "DMSO3_lane1.count",   
                               "Aph_KH7_1.count", "Aph_KH7_2.count", "Aph_KH7_3.count" ))
dim(genes.counts.Aph_KH7)
head(genes.counts.Aph_KH7)

###############################################################################
###############################################################################

genes.counts.KH7 <- subset(genes.counts, select=c( 
                           "DMSO1_lane1.count", 
                           "DMSO2_lane1.count", 
                           "DMSO3_lane1.count",  
                           "KH7_1.count", "KH7_2.count", "KH7_3.count" ))
dim(genes.counts.KH7)
head(genes.counts.KH7)

###############################################################################
###############################################################################

genes.counts.Noc <- subset(genes.counts, select=c( 
                           "DMSO1_lane1.count", 
                           "DMSO2_lane1.count", 
                           "DMSO3_lane1.count",  
                           "Noc_1.count", "Noc_2.count", "Noc_3.count" ))
dim(genes.counts.Noc)
head(genes.counts.Noc)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

############################# STARTING TO ASSESS THE DIFFERENTIAL EXPRESSION : using LIMMA :

############################################################################################
############################################################################################
############################################################################################
############################################################################################
################################################## using LIMMA for each individual MATRIX : 

#### genes.counts.Aph 
#### genes.counts.Aph_KH7
#### genes.counts.KH7
#### genes.counts.Noc

############################################################################################
############################################################################################
############################################################################################
############################################################################################

#### having a model in an OLD PIECE of CODE :

#### setting up the groups and the subjects
#### group <- factor(c("csc","csc","csc","csc","csc","non","non","non","non","non"))
#### subject <- factor(c(1,2,3,4,5,1,2,3,4,5))

#### setting up the design and the contrast matrix
#### design <- model.matrix(~0+group+subject)
#### contrast.matrix <- makeContrasts(groupcsc-groupnon, levels=design)

############################################################################################
############################################################################################
############################################################################################
############################################################################################

eset <- genes.counts.Aph
eset_name <- deparse(substitute(genes.counts.Aph)) ### in order to get the name of the DF

#### genes.counts.Aph
 
group <- factor(c("DMSO", "DMSO", "DMSO", "Aph", "Aph","Aph"))
subject <- factor(c(1,2,3, 1,2,3))

#### setting up the design and the contrast matrix

design <- model.matrix(~0+group+subject)
contrast.matrix <- makeContrasts(groupAph-groupDMSO, levels=design)

design
contrast.matrix

############################################################################################
############################################################################################

### filtering the genes based on CPM :

y <- DGEList(counts=eset, group=group)

### keep <- rowSums(cpm(y, lib.size=libsize)>1) >= 3

keep <- rowSums( cpm(y) > 0.5) >= 6
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

####################################################### we can use y$counts for PCA analysis

### computing the normalization factors :
y <- calcNormFactors(y)

### using the VOOM transformation :
v <- voom(y, design, plot=FALSE)

### the LINEAR FIT in LIMMA :
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

### obtaining and writing the results :
results_limma <- topTable(fit2, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns 

results_limma$Gene <- rownames(results_limma)

### separating the names of the GENES into 1st_PART and NUMBER :

results_limma$GENE <- results_limma$Gene 

results_limma.sep <- separate(data=results_limma, col=Gene, into = c("Gene", "ID"), sep = ":")

head(results_limma.sep)
dim(results_limma.sep)

### writing the results to a file :

write.table(results_limma.sep, file=paste("analysis.LIMMA.", eset_name, sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 :

results_limma.deg <- results_limma.sep[results_limma.sep$adj.P.Val < 0.05,]

head(results_limma.deg)
dim(results_limma.deg)

write.table(results_limma.deg, file=paste("analysis.LIMMA.", eset_name, ".only.DEG", sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)


### computing the number of DEG for FDR < 0.05 and FC > 1.2  : UP-REGULATED GENES :

results_limma.deg.up <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                          (results_limma.sep$logFC > log2(1.2) )  ,]

head(results_limma.deg.up)
dim(results_limma.deg.up)

write.table(results_limma.deg.up, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP", sep=""),
                                  sep="\t", 
                                  quote=FALSE, eol="\n", 
                                  row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC < -1.2 : DOWN-REGULATED GENES :

results_limma.deg.down <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                            (results_limma.sep$logFC < -log2(1.2) )  ,]

head(results_limma.deg.down)
dim(results_limma.deg.down)

write.table(results_limma.deg.down, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN", sep=""),
                                    sep="\t", 
                                    quote=FALSE, eol="\n", 
                                    row.names=FALSE, col.names=TRUE)


############################################################################################
############################################################################################
############################################### saving the results into another DATAFRAME  :

eset <- genes.counts.Aph
eset_name <- deparse(substitute(genes.counts.Aph)) ##### in order to get the name of the DF
genes.counts.Aph.results.limma <- results_limma.sep

head(genes.counts.Aph.results.limma)
dim(genes.counts.Aph.results.limma)

genes.counts.Aph.results.limma.deg <- results_limma.deg
genes.counts.Aph.results.limma.deg.up <- results_limma.deg.up
genes.counts.Aph.results.limma.deg.down <- results_limma.deg.down

dim(genes.counts.Aph.results.limma.deg.up)
dim(genes.counts.Aph.results.limma.deg.down)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

eset <- genes.counts.Aph_KH7
eset_name <- deparse(substitute(genes.counts.Aph_KH7))

#### genes.counts.Aph_KH7
 
group <- factor(c("DMSO", "DMSO", "DMSO", "Aph_KH7", "Aph_KH7", "Aph_KH7"))
subject <- factor(c(1,2,3, 1,2,3))

### setting up the design and the contrast matrix

design <- model.matrix(~0+group+subject)
contrast.matrix <- makeContrasts(groupAph_KH7-groupDMSO, levels=design)

design
contrast.matrix

############################################################################################

### filtering the genes based on CPM :
y <- DGEList(counts=eset, group=group)

### keep <- rowSums(cpm(y, lib.size=libsize)>1) >= 3
keep <- rowSums( cpm(y) > 0.5) >= 6
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

####################################################### we can use y$counts for PCA analysis

### computing the normalization factors :
y <- calcNormFactors(y)

### using the VOOM transformation :
v <- voom(y, design, plot=FALSE)

### the LINEAR FIT in LIMMA :
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

### obtaining and writing the results :
results_limma <- topTable(fit2, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns 

results_limma$Gene <- rownames(results_limma)

### separating the names of the GENES into 1st_PART and NUMBER :

results_limma$GENE <- results_limma$Gene 

results_limma.sep <- separate(data=results_limma, col=Gene, into = c("Gene", "ID"), sep = ":")

head(results_limma.sep)
dim(results_limma.sep)

### writing the results to a file :

write.table(results_limma.sep, file=paste("analysis.LIMMA.", eset_name, sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 :

results_limma.deg <- results_limma.sep[results_limma.sep$adj.P.Val < 0.05,]

head(results_limma.deg)
dim(results_limma.deg)

write.table(results_limma.deg, file=paste("analysis.LIMMA.", eset_name, ".only.DEG", sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC > 1.2  : UP-REGULATED GENES :

results_limma.deg.up <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                          (results_limma.sep$logFC > log2(1.2) )  ,]

head(results_limma.deg.up)
dim(results_limma.deg.up)

write.table(results_limma.deg.up, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP", sep=""),
                                  sep="\t", 
                                  quote=FALSE, eol="\n", 
                                  row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC < -1.2 : DOWN-REGULATED GENES :

results_limma.deg.down <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                            (results_limma.sep$logFC < -log2(1.2) )  ,]

head(results_limma.deg.down)
dim(results_limma.deg.down)

write.table(results_limma.deg.down, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN", sep=""),
                                    sep="\t", 
                                    quote=FALSE, eol="\n", 
                                    row.names=FALSE, col.names=TRUE)

############################################################################################
############################################################################################
### saving the results into another DATAFRAME to be used at a later time point :

eset <- genes.counts.Aph_KH7
eset_name <- deparse(substitute(genes.counts.Aph_KH7))
genes.counts.Aph_KH7.results.limma <- results_limma.sep

head(genes.counts.Aph_KH7.results.limma)
dim(genes.counts.Aph_KH7.results.limma)

genes.counts.Aph_KH7.results.limma.deg <- results_limma.deg
genes.counts.Aph_KH7.results.limma.deg.up <- results_limma.deg.up
genes.counts.Aph_KH7.results.limma.deg.down <- results_limma.deg.down 

dim(genes.counts.Aph_KH7.results.limma.deg.up)
dim(genes.counts.Aph_KH7.results.limma.deg.down)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

eset <- genes.counts.KH7
eset_name <- deparse(substitute(genes.counts.KH7))

#### genes.counts.KH7
 
group <- factor(c("DMSO","DMSO","DMSO", "KH7","KH7","KH7"))
subject <- factor(c(1,2,3, 1,2,3))

### setting up the design and the contrast matrix

design <- model.matrix(~0+group+subject)
contrast.matrix <- makeContrasts(groupKH7-groupDMSO, levels=design)

design
contrast.matrix

############################################################################################

### filtering the genes based on CPM :
y <- DGEList(counts=eset, group=group)

### keep <- rowSums(cpm(y, lib.size=libsize)>1) >= 3
keep <- rowSums( cpm(y) > 0.5) >= 6
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

####################################################### we can use y$counts for PCA analysis

### computing the normalization factors :
y <- calcNormFactors(y)

### using the VOOM transformation :
v <- voom(y, design, plot=FALSE)

### the LINEAR FIT in LIMMA :
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

### obtaining and writing the results :
results_limma <- topTable(fit2, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns 

results_limma$Gene <- rownames(results_limma)

### separating the names of the GENES into 1st_PART and NUMBER :

results_limma$GENE <- results_limma$Gene 

results_limma.sep <- separate(data=results_limma, col=Gene, into = c("Gene", "ID"), sep = ":")

head(results_limma.sep)
dim(results_limma.sep)

### writing the results to a file :

write.table(results_limma.sep, file=paste("analysis.LIMMA.", eset_name, sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 :

results_limma.deg <- results_limma.sep[results_limma.sep$adj.P.Val < 0.05,]

head(results_limma.deg)
dim(results_limma.deg)

write.table(results_limma.deg, file=paste("analysis.LIMMA.", eset_name, ".only.DEG", sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC > 1.2  : UP-REGULATED GENES :

results_limma.deg.up <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                          (results_limma.sep$logFC > log2(1.2) )  ,]

head(results_limma.deg.up)
dim(results_limma.deg.up)

write.table(results_limma.deg.up, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP", sep=""),
                                  sep="\t", 
                                  quote=FALSE, eol="\n", 
                                  row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC < -1.2 : DOWN-REGULATED GENES :

results_limma.deg.down <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                            (results_limma.sep$logFC < -log2(1.2) )  ,]

head(results_limma.deg.down)
dim(results_limma.deg.down)

write.table(results_limma.deg.down, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN", sep=""),
                                    sep="\t", 
                                    quote=FALSE, eol="\n", 
                                    row.names=FALSE, col.names=TRUE)

############################################################################################
############################################################################################

eset <- genes.counts.KH7
eset_name <- deparse(substitute(genes.counts.KH7))
genes.counts.KH7.results.limma <- results_limma.sep

head(genes.counts.KH7.results.limma)
dim(genes.counts.KH7.results.limma)

genes.counts.KH7.results.limma.deg <- results_limma.deg
genes.counts.KH7.results.limma.deg.up <- results_limma.deg.up
genes.counts.KH7.results.limma.deg.down <- results_limma.deg.down 

dim(genes.counts.KH7.results.limma.deg.up)
dim(genes.counts.KH7.results.limma.deg.down)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

eset <- genes.counts.Noc
eset_name <- deparse(substitute(genes.counts.Noc))

#### genes.counts.Noc
 
group <- factor(c("DMSO", "DMSO", "DMSO", "Noc", "Noc", "Noc"))
subject <- factor(c(1,2,3, 1,2,3))

### setting up the design and the contrast matrix

design <- model.matrix(~0+group+subject)
contrast.matrix <- makeContrasts(groupNoc-groupDMSO, levels=design)

design
contrast.matrix

############################################################################################

### filtering the genes based on CPM :
y <- DGEList(counts=eset, group=group)

### keep <- rowSums(cpm(y, lib.size=libsize)>1) >= 3
keep <- rowSums( cpm(y) > 0.5) >= 6
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

####################################################### we can use y$counts for PCA analysis

### computing the normalization factors :
y <- calcNormFactors(y)

### using the VOOM transformation :
v <- voom(y, design, plot=FALSE)

### the LINEAR FIT in LIMMA :
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

### obtaining and writing the results :
results_limma <- topTable(fit2, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns 

results_limma$Gene <- rownames(results_limma)

### separating the names of the GENES into 1st_PART and NUMBER :

results_limma$GENE <- results_limma$Gene 

results_limma.sep <- separate(data=results_limma, col=Gene, into = c("Gene", "ID"), sep = ":")

head(results_limma.sep)
dim(results_limma.sep)

### writing the results to a file :

write.table(results_limma.sep, file=paste("analysis.LIMMA.", eset_name, sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 :

results_limma.deg <- results_limma.sep[results_limma.sep$adj.P.Val < 0.05,]

head(results_limma.deg)
dim(results_limma.deg)

write.table(results_limma.deg, file=paste("analysis.LIMMA.", eset_name, ".only.DEG", sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC > 1.2  : UP-REGULATED GENES :

results_limma.deg.up <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                          (results_limma.sep$logFC > log2(1.2) )  ,]

head(results_limma.deg.up)
dim(results_limma.deg.up)

write.table(results_limma.deg.up, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP", sep=""),
                                  sep="\t", 
                                  quote=FALSE, eol="\n", 
                                  row.names=FALSE, col.names=TRUE)

### computing the number of DEG for FDR < 0.05 and FC < -1.2 : DOWN-REGULATED GENES :

results_limma.deg.down <- results_limma.sep[(results_limma.sep$adj.P.Val < 0.05) & 
                                            (results_limma.sep$logFC < -log2(1.2) )  ,]

head(results_limma.deg.down)
dim(results_limma.deg.down)

write.table(results_limma.deg.down, file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN", sep=""),
                                    sep="\t", 
                                    quote=FALSE, eol="\n", 
                                    row.names=FALSE, col.names=TRUE)

############################################################################################

eset <- genes.counts.Noc
eset_name <- deparse(substitute(genes.counts.Noc))
genes.counts.Noc.results.limma <- results_limma.sep

head(genes.counts.Noc.results.limma)
dim(genes.counts.Noc.results.limma)

genes.counts.Noc.results.limma.deg <- results_limma.deg
genes.counts.Noc.results.limma.deg.up <- results_limma.deg.up
genes.counts.Noc.results.limma.deg.down <- results_limma.deg.down 

dim(genes.counts.Noc.results.limma.deg.up)
dim(genes.counts.Noc.results.limma.deg.down)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

#### AT THIS MOMENT, we would like to INTEGRATE all the DATAFILES from LIMMA that we have :

#### genes OR genes.counts 

#### genes.counts.Aph.results.limma 
#### genes.counts.Aph_KH7.results.limma
#### genes.counts.KH7.results.limma
#### genes.counts.Noc.results.limma

dim(genes)
dim(genes.counts.Aph.results.limma)
dim(genes.counts.Aph_KH7.results.limma)
dim(genes.counts.KH7.results.limma)
dim(genes.counts.Noc.results.limma)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

################################## we will have to change the names of columns, 
################################## because all the dataframes have the same COLUMN NAMES

################################## working with a DATAFRAME : Aph 

colnames(genes.counts.Aph.results.limma)[1] <- paste("logFC" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[2] <- paste("AveExpr" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[3] <- paste("t" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[4] <- paste("P.Value" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[5] <- paste("adj.P.Val" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[6] <- paste("B" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[7] <- paste("Gene" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)[8] <- paste("ID" ,"Aph" , sep=":")
colnames(genes.counts.Aph.results.limma)
# colnames(genes.counts.Aph.results.limma)[9] <- paste("GENE" ,"Aph" , sep=":")

################################## working with a DATAFRAME : KH7

colnames(genes.counts.KH7.results.limma)[1] <- paste("logFC" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[2] <- paste("AveExpr" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[3] <- paste("t" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[4] <- paste("P.Value" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[5] <- paste("adj.P.Val" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[6] <- paste("B" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[7] <- paste("Gene" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)[8] <- paste("ID" ,"KH7" , sep=":")
colnames(genes.counts.KH7.results.limma)
# colnames(genes.counts.KH7.results.limma)[9] <- paste("GENE" ,"KH7" , sep=":")

################################## working with a DATAFRAME : Aph_KH7

colnames(genes.counts.Aph_KH7.results.limma)[1] <- paste("logFC" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[2] <- paste("AveExpr" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[3] <- paste("t" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[4] <- paste("P.Value" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[5] <- paste("adj.P.Val" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[6] <- paste("B" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[7] <- paste("Gene" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)[8] <- paste("ID" ,"Aph_KH7" , sep=":")
colnames(genes.counts.Aph_KH7.results.limma)
# colnames(genes.counts.Aph_KH7.results.limma)[9] <- paste("GENE" ,"Aph_KH7" , sep=":")

################################## working with a DATAFRAME : Noc

colnames(genes.counts.Noc.results.limma)[1] <- paste("logFC" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[2] <- paste("AveExpr" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[3] <- paste("t" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[4] <- paste("P.Value" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[5] <- paste("adj.P.Val" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[6] <- paste("B" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[7] <- paste("Gene" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)[8] <- paste("ID" ,"Noc" , sep=":")
colnames(genes.counts.Noc.results.limma)
# colnames(genes.counts.Noc.results.limma)[9] <- paste("GENE" ,"Noc" , sep=":")

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

### now integrating these data structures ; we can make DATA TABLES :

genes.dt <- as.data.table(genes)

genes.counts.Aph.results.limma.dt <- as.data.table(genes.counts.Aph.results.limma)
genes.counts.Aph_KH7.results.limma.dt <- as.data.table(genes.counts.Aph_KH7.results.limma)
genes.counts.KH7.results.limma.dt <- as.data.table(genes.counts.KH7.results.limma)
genes.counts.Noc.results.limma.dt <- as.data.table(genes.counts.Noc.results.limma)

#################################### setting up the KEYS :

setkeyv(genes.dt, c('GENE_NAME_ID'))

setkeyv(genes.counts.Aph.results.limma.dt, c('GENE'))
setkeyv(genes.counts.Aph_KH7.results.limma.dt, c('GENE'))
setkeyv(genes.counts.KH7.results.limma.dt, c('GENE'))
setkeyv(genes.counts.Noc.results.limma.dt, c('GENE'))

#########################################################################################################
#########################################################################################################
#########################################################################################################

integration.all.samples.dt <- genes.dt[genes.counts.Aph.results.limma.dt,][genes.counts.Aph_KH7.results.limma.dt,][genes.counts.KH7.results.limma.dt,][genes.counts.Noc.results.limma.dt,]

head(integration.all.samples.dt)
dim(integration.all.samples.dt)

write.table(integration.all.samples.dt, file=paste("analysis.LIMMA.integrating.all.samples.with.data.table", sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

################### starting from the dataframe "genes" :

dim(genes)
head(genes)

################### integrate with : genes.counts.Aph.results.limma

integration.step1 <- merge(genes, 
                           genes.counts.Aph.results.limma, 
                           by.x = "GENE_NAME_ID" , 
                           by.y = "GENE",  
                           all.x = TRUE)

head(integration.step1)
dim(integration.step1)

################### integrate with : genes.counts.Aph_KH7.results.limma

integration.step2 <- merge(integration.step1, 
                           genes.counts.Aph_KH7.results.limma, 
                           by.x = "GENE_NAME_ID" , 
                           by.y = "GENE",  
                           all.x = TRUE)

head(integration.step2)
dim(integration.step2)

################### integrate with : genes.counts.KH7.results.limma

integration.step3 <- merge(integration.step2, 
                           genes.counts.KH7.results.limma, 
                           by.x = "GENE_NAME_ID" , 
                           by.y = "GENE",  
                           all.x = TRUE)

head(integration.step3)
dim(integration.step3)

################### integrate with : genes.counts.Noc.results.limma

integration.step4 <- merge(integration.step3, 
                           genes.counts.Noc.results.limma, 
                           by.x = "GENE_NAME_ID" , 
                           by.y = "GENE",  
                           all.x = TRUE)

head(integration.step4)
dim(integration.step4)

############################################################################################
############################################################################################
########################## we are going to write the file for computing the fold changes ...

write.table(integration.step4, file=paste("analysis.LIMMA.integrating.all.samples.all.genes.in.4.STEPS", sep=""),
                               sep="\t", 
                               quote=FALSE, eol="\n", 
                               row.names=FALSE, col.names=TRUE)

############################################################################################
############################################################################################
 
########### starting from the dataframe "integration.step4", to separate the DEG, function of FDR, and FC :

# x <- integration.step4 

# colnames(x)
# [1] "GENE_NAME_ID"      "CHR"               "START"            
# [4] "END"               "STRAND"            "GENE_ID"          
# [7] "GENE_NAME"         "GENE_TYPE"         "DMSO1_lane1.count"
#[10] "DMSO1_lane1.TPM"   "DMSO1_lane1.FPKM"  "DMSO1_lane2.count"
#[13] "DMSO1_lane2.TPM"   "DMSO1_lane2.FPKM"  "DMSO2_lane1.count"
#[16] "DMSO2_lane1.TPM"   "DMSO2_lane1.FPKM"  "DMSO2_lane2.count"
#[19] "DMSO2_lane2.TPM"   "DMSO2_lane2.FPKM"  "DMSO3_lane1.count"
#[22] "DMSO3_lane1.TPM"   "DMSO3_lane1.FPKM"  "DMSO3_lane2.count"
#[25] "DMSO3_lane2.TPM"   "DMSO3_lane2.FPKM"  "Aph1.count"       
#[28] "Aph1.TPM"          "Aph1.FPKM"         "Aph2.count"       
#[31] "Aph2.TPM"          "Aph2.FPKM"         "Aph3.count"       
#[34] "Aph3.TPM"          "Aph3.FPKM"         "Aph_KH7_1.count"  
#[37] "Aph_KH7_1.TPM"     "Aph_KH7_1.FPKM"    "Aph_KH7_2.count"  
#[40] "Aph_KH7_2.TPM"     "Aph_KH7_2.FPKM"    "Aph_KH7_3.count"  
#[43] "Aph_KH7_3.TPM"     "Aph_KH7_3.FPKM"    "KH7_1.count"      
#[46] "KH7_1.TPM"         "KH7_1.FPKM"        "KH7_2.count"      
#[49] "KH7_2.TPM"         "KH7_2.FPKM"        "KH7_3.count"      
#[52] "KH7_3.TPM"         "KH7_3.FPKM"        "Noc_1.count"      
#[55] "Noc_1.TPM"         "Noc_1.FPKM"        "Noc_2.count"      
#[58] "Noc_2.TPM"         "Noc_2.FPKM"        "Noc_3.count"      
#[61] "Noc_3.TPM"         "Noc_3.FPKM"        "ID"               
#[64] "logFC:Aph"         "AveExpr:Aph"       "t:Aph"            
#[67] "P.Value:Aph"       "adj.P.Val:Aph"     "B:Aph"            
#[70] "Gene:Aph"          "ID:Aph"            "logFC:Aph_KH7"    
#[73] "AveExpr:Aph_KH7"   "t:Aph_KH7"         "P.Value:Aph_KH7"  
#[76] "adj.P.Val:Aph_KH7" "B:Aph_KH7"         "Gene:Aph_KH7"     
#[79] "ID:Aph_KH7"        "logFC:KH7"         "AveExpr:KH7"      
#[82] "t:KH7"             "P.Value:KH7"       "adj.P.Val:KH7"    
#[85] "B:KH7"             "Gene:KH7"          "ID:KH7"           
#[88] "logFC:Noc"         "AveExpr:Noc"       "t:Noc"            
#[91] "P.Value:Noc"       "adj.P.Val:Noc"     "B:Noc"            
#[94] "Gene:Noc"          "ID:Noc"  

#############################################################################################
#############################################################################################
############### using as criteria for filtering the following fields :

# "DMSO1_lane1.FPKM"  
# "DMSO1_lane2.FPKM"  
# "DMSO2_lane1.FPKM"  
# "DMSO2_lane2.FPKM"  
# "DMSO3_lane1.FPKM"  
# "DMSO3_lane2.FPKM"         
# "Aph1.FPKM"                
# "Aph2.FPKM"                
# "Aph3.FPKM"           
# "Aph_KH7_1.FPKM"      
# "Aph_KH7_2.FPKM"      
# "Aph_KH7_3.FPKM"         
# "KH7_1.FPKM"              
# "KH7_2.FPKM"              
# "KH7_3.FPKM"              
# "Noc_1.FPKM"              
# "Noc_2.FPKM"

#############################################################################################
################ considering the FPKM, FC, and FDR :
#############################################################################################

#"logFC:Aph"                     
#"adj.P.Val:Aph"       
#"logFC:Aph_KH7"    
#"adj.P.Val:Aph_KH7"    
#"logFC:KH7"
#"adj.P.Val:KH7"
#"logFC:Noc"
#"adj.P.Val:Noc"

x <- integration.step4

head(x$"logFC:Aph")                     
head(x$"adj.P.Val:Aph")       
head(x$"logFC:Aph_KH7")    
head(x$"adj.P.Val:Aph_KH7")    
head(x$"logFC:KH7")
head(x$"adj.P.Val:KH7")
head(x$"logFC:Noc")
head(x$"adj.P.Val:Noc")
         
head(x)
tail(x)
dim(x)
 
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

### to integrate these DATAFRAMES with X, and to print it :

############################################################################################
############################################################################################
####################### considering the comparisons DMSO:Aph :

eset_name <- deparse(substitute(genes.counts.Aph))

# genes.counts.Aph.results.limma.deg 
# genes.counts.Aph.results.limma.deg.up 
# genes.counts.Aph.results.limma.deg.down 

genes.counts.Aph.results.limma.deg.up.and.x <- merge(genes.counts.Aph.results.limma.deg.up, 
                                                     x,
                                                     by.x="GENE",
                                                     by.y="GENE_NAME_ID",
                                                     all.x = TRUE) 


write.table(genes.counts.Aph.results.limma.deg.up.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)
     
head(genes.counts.Aph.results.limma.deg.up.and.x)
tail(genes.counts.Aph.results.limma.deg.up.and.x)
dim(genes.counts.Aph.results.limma.deg.up.and.x) 

genes.counts.Aph.results.limma.deg.down.and.x <- merge(genes.counts.Aph.results.limma.deg.down, 
                                                     x,
                                                     by.x="GENE",
                                                     by.y="GENE_NAME_ID",
                                                     all.x = TRUE) 

write.table(genes.counts.Aph.results.limma.deg.down.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)

head(genes.counts.Aph.results.limma.deg.down.and.x)
tail(genes.counts.Aph.results.limma.deg.down.and.x)
dim(genes.counts.Aph.results.limma.deg.down.and.x) 

############################################################################################
############################################################################################
####################### considering the comparisons DMSO:Aph_KH7 :

eset_name <- deparse(substitute(genes.counts.Aph_KH7))

# genes.counts.Aph_KH7.results.limma.deg 
# genes.counts.Aph_KH7.results.limma.deg.up 
# genes.counts.Aph_KH7.results.limma.deg.down 

genes.counts.Aph_KH7.results.limma.deg.up.and.x <- merge(genes.counts.Aph_KH7.results.limma.deg.up, 
                                                         x,
                                                         by.x="GENE",
                                                         by.y="GENE_NAME_ID",
                                                         all.x = TRUE) 

write.table(genes.counts.Aph_KH7.results.limma.deg.up.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)

head(genes.counts.Aph_KH7.results.limma.deg.up.and.x)
tail(genes.counts.Aph_KH7.results.limma.deg.up.and.x)
dim(genes.counts.Aph_KH7.results.limma.deg.up.and.x)

genes.counts.Aph_KH7.results.limma.deg.down.and.x <- merge(genes.counts.Aph_KH7.results.limma.deg.down, 
                                                         x,
                                                         by.x="GENE",
                                                         by.y="GENE_NAME_ID",
                                                         all.x = TRUE)

write.table(genes.counts.Aph_KH7.results.limma.deg.down.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)

head(genes.counts.Aph_KH7.results.limma.deg.down.and.x)
tail(genes.counts.Aph_KH7.results.limma.deg.down.and.x)
dim(genes.counts.Aph_KH7.results.limma.deg.down.and.x)

############################################################################################
############################################################################################
####################### considering the comparisons DMSO:KH7 :

eset_name <- deparse(substitute(genes.counts.KH7))

# genes.counts.KH7.results.limma.deg 
# genes.counts.KH7.results.limma.deg.up 
# genes.counts.KH7.results.limma.deg.down 

genes.counts.KH7.results.limma.deg.up.and.x <- merge(genes.counts.KH7.results.limma.deg.up, 
                                                         x,
                                                         by.x="GENE",
                                                         by.y="GENE_NAME_ID",
                                                         all.x = TRUE)

write.table(genes.counts.KH7.results.limma.deg.up.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE) 

head(genes.counts.KH7.results.limma.deg.up.and.x)
tail(genes.counts.KH7.results.limma.deg.up.and.x)
dim(genes.counts.KH7.results.limma.deg.up.and.x)

genes.counts.KH7.results.limma.deg.down.and.x  <- merge(genes.counts.KH7.results.limma.deg.down, 
                                                         x,
                                                         by.x="GENE",
                                                         by.y="GENE_NAME_ID",
                                                         all.x = TRUE)

write.table(genes.counts.KH7.results.limma.deg.down.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)

head(genes.counts.KH7.results.limma.deg.down.and.x)
tail(genes.counts.KH7.results.limma.deg.down.and.x)
dim(genes.counts.KH7.results.limma.deg.down.and.x)

############################################################################################
############################################################################################
####################### considering the comparisons DMSO:Noc :

eset_name <- deparse(substitute(genes.counts.Noc))

# genes.counts.Noc.results.limma.deg 
# genes.counts.Noc.results.limma.deg.up 
# genes.counts.Noc.results.limma.deg.down 


genes.counts.Noc.results.limma.deg.up.and.x <- merge(genes.counts.Noc.results.limma.deg.up, 
                                                         x,
                                                         by.x="GENE",
                                                         by.y="GENE_NAME_ID",
                                                         all.x = TRUE) 

write.table(genes.counts.Noc.results.limma.deg.up.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.UP.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)

head(genes.counts.Noc.results.limma.deg.up.and.x)
tail(genes.counts.Noc.results.limma.deg.up.and.x)
dim(genes.counts.Noc.results.limma.deg.up.and.x)

genes.counts.Noc.results.limma.deg.down.and.x  <- merge(genes.counts.Noc.results.limma.deg.down, 
                                                         x,
                                                         by.x="GENE",
                                                         by.y="GENE_NAME_ID",
                                                         all.x = TRUE)

write.table(genes.counts.Noc.results.limma.deg.down.and.x, 
            file=paste("analysis.LIMMA.", eset_name, ".only.DEG.and.DOWN.info.all.samples", sep=""),
            sep="\t", 
            quote=FALSE, eol="\n", 
            row.names=FALSE, col.names=TRUE)

head(genes.counts.Noc.results.limma.deg.down.and.x)
tail(genes.counts.Noc.results.limma.deg.down.and.x)
dim(genes.counts.Noc.results.limma.deg.down.and.x)

############################################################################################
############################################################################################
############################################################################################
