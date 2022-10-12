############################################################################################
############################################################################################
############################################################################################
############################################################################################

library("ggplot2")
library("reshape2")
library("data.table")

library("limma")
library("Glimma")
library("edgeR")
library("DESeq2")

library("pheatmap")
library("ComplexHeatmap")
library("gplots")
library("scatterplot3d")
library("enrichR")
library("tidyr")
library("plyr")
library("dplyr")
library("RColorBrewer")

library("VennDiagram") 
library("Vennerable")
library("gplots")
library("data.table")

############################################################################################
############################################################################################
############################################################################################
############################################################################################

# starting with the OUTPUT files from RSEM

SAMPLE1 = "Control_1" 
SAMPLE2 = "Control_2" 
SAMPLE3 = "Sox11_1"
SAMPLE4 = "Sox11_2"
SAMPLE5 = "Sox11MUT_1"
SAMPLE6 = "Sox11MUT_2"

# the files that we are reading are the following :

# gencode.vM4.annotation-tRNAs-ERCC.gtf.processed.unique.gene.only.ENSMUSG.txt
# control1_anno_rsem.genes.results
# control2_anno_rsem.genes.results
# SOX11_1_anno_rsem.genes.results
# SOX11_2_anno_rsem.genes.results
# SOX11_mutant1_anno_rsem.genes.results
# SOX11_mutant2_anno_rsem.genes.results


############################################################################################
############################################################################################
############################################################################################
############################################################################################

genes <- read.delim("gencode.vM4.annotation-tRNAs-ERCC.gtf.processed.unique.gene.only.ENSMUSG.txt", 
                     sep="\t", header=T, stringsAsFactors=F)

head(genes)
dim(genes)

genes.dt <- as.data.table(genes)
head(genes.dt)
dim(genes.dt)  

###### to integrate these files : reading the files and changing the names of the columns

name <- "gencode.vM4.annotation-tRNAs-ERCC.gtf.processed.unique.gene.only.ENSMUSG.txt"

############################################################################################
############################################################################################

# The file with the genes list has the following FORMAT :

# [1] "gene_id"           "gene_name"         "gene_type"        
# [4] "gene_status"       "transcript_id"     "transcript_type"  
# [7] "transcript_status" "transcript_name"  

# The files with the data from RSEM have the following FORMAT :

# [1] "gene_id"                               
# [2] "transcript_id.s."                      
# [3] "length"                                
# [4] "effective_length"                      
# [5] "expected_count"                        
# [6] "TPM"                                   
# [7] "FPKM"                                  
# [8] "posterior_mean_count"                  
# [9] "posterior_standard_deviation_of_count" 
# [10] "pme_TPM"                               
# [11] "pme_FPKM"                              
# [12] "TPM_ci_lower_bound"                    
# [13] "TPM_ci_upper_bound"                    
# [14] "TPM_coefficient_of_quartile_variation" 
# [15] "FPKM_ci_lower_bound"                   
# [16] "FPKM_ci_upper_bound"                   
# [17] "FPKM_coefficient_of_quartile_variation"

############################################################################################
############################################################################################
############################################################################################ SAMPLE1

SAMPLE1.file <- read.delim("control1_anno_rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

SAMPLE1.df   <- data.frame(matrix(ncol = 4, nrow = dim(SAMPLE1.file)[1]))

colnames(SAMPLE1.df) <- c(paste(SAMPLE1, "gene", sep="."), 
                          paste(SAMPLE1, "count", sep="."), 
                          paste(SAMPLE1, "TPM", sep="."), 
                          paste(SAMPLE1, "FPKM", sep="."))

SAMPLE1.df[,1] = SAMPLE1.file$gene_id
SAMPLE1.df[,2] = SAMPLE1.file$expected_count
SAMPLE1.df[,3] = SAMPLE1.file$TPM
SAMPLE1.df[,4] = SAMPLE1.file$FPKM
                                 
dim(SAMPLE1.df)

# head(SAMPLE1.df)
# write.table(SAMPLE1.df, file=paste(SAMPLE1, "output.txt", sep="."), 
#                         quote = FALSE, sep = "\t", 
#                         row.names = FALSE, col.names = TRUE)
             
############################################################################################
############################################################################################
############################################################################################ SAMPLE2

SAMPLE2.file <- read.delim("control2_anno_rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

SAMPLE2.df   <- data.frame(matrix(ncol = 4, nrow = dim(SAMPLE2.file)[1]))

colnames(SAMPLE2.df) <- c(paste(SAMPLE2, "gene", sep="."), 
                          paste(SAMPLE2, "count", sep="."), 
                          paste(SAMPLE2, "TPM", sep="."), 
                          paste(SAMPLE2, "FPKM", sep="."))

SAMPLE2.df[,1] = SAMPLE2.file$gene_id
SAMPLE2.df[,2] = SAMPLE2.file$expected_count
SAMPLE2.df[,3] = SAMPLE2.file$TPM
SAMPLE2.df[,4] = SAMPLE2.file$FPKM
                                 
dim(SAMPLE2.df)

# head(SAMPLE2.df)
# write.table(SAMPLE2.df, file=paste(SAMPLE2, "output.txt", sep="."), 
#                         quote = FALSE, sep = "\t", 
#                         row.names = FALSE, col.names = TRUE)

############################################################################################
############################################################################################
############################################################################################ SAMPLE3

SAMPLE3.file <- read.delim("SOX11_1_anno_rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

SAMPLE3.df   <- data.frame(matrix(ncol = 4, nrow = dim(SAMPLE3.file)[1]))

colnames(SAMPLE3.df) <- c(paste(SAMPLE3, "gene", sep="."), 
                          paste(SAMPLE3, "count", sep="."), 
                          paste(SAMPLE3, "TPM", sep="."), 
                          paste(SAMPLE3, "FPKM", sep="."))

SAMPLE3.df[,1] = SAMPLE3.file$gene_id
SAMPLE3.df[,2] = SAMPLE3.file$expected_count
SAMPLE3.df[,3] = SAMPLE3.file$TPM
SAMPLE3.df[,4] = SAMPLE3.file$FPKM
                                 
dim(SAMPLE3.df)

# head(SAMPLE3.df)
# write.table(SAMPLE3.df, file=paste(SAMPLE3, "output.txt", sep="."), 
#                         quote = FALSE, sep = "\t", 
#                         row.names = FALSE, col.names = TRUE)

############################################################################################
############################################################################################
############################################################################################ SAMPLE4

SAMPLE4.file <- read.delim("SOX11_2_anno_rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

SAMPLE4.df   <- data.frame(matrix(ncol = 4, nrow = dim(SAMPLE4.file)[1]))

colnames(SAMPLE4.df) <- c(paste(SAMPLE4, "gene", sep="."), 
                          paste(SAMPLE4, "count", sep="."), 
                          paste(SAMPLE4, "TPM", sep="."), 
                          paste(SAMPLE4, "FPKM", sep="."))

SAMPLE4.df[,1] = SAMPLE4.file$gene_id
SAMPLE4.df[,2] = SAMPLE4.file$expected_count
SAMPLE4.df[,3] = SAMPLE4.file$TPM
SAMPLE4.df[,4] = SAMPLE4.file$FPKM
                                 
dim(SAMPLE4.df)

# head(SAMPLE4.df)
# write.table(SAMPLE4.df, file=paste(SAMPLE4, "output.txt", sep="."), 
#                         quote = FALSE, sep = "\t", 
#                         row.names = FALSE, col.names = TRUE)

############################################################################################
############################################################################################
############################################################################################ SAMPLE5

SAMPLE5.file <- read.delim("SOX11_mutant1_anno_rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

SAMPLE5.df   <- data.frame(matrix(ncol = 4, nrow = dim(SAMPLE5.file)[1]))

colnames(SAMPLE5.df) <- c(paste(SAMPLE5, "gene", sep="."), 
                          paste(SAMPLE5, "count", sep="."), 
                          paste(SAMPLE5, "TPM", sep="."), 
                          paste(SAMPLE5, "FPKM", sep="."))

SAMPLE5.df[,1] = SAMPLE5.file$gene_id
SAMPLE5.df[,2] = SAMPLE5.file$expected_count
SAMPLE5.df[,3] = SAMPLE5.file$TPM
SAMPLE5.df[,4] = SAMPLE5.file$FPKM
                                 
dim(SAMPLE5.df)

# head(SAMPLE5.df)
# write.table(SAMPLE5.df, file=paste(SAMPLE5, "output.txt", sep="."), 
#                         quote = FALSE, sep = "\t", 
#                         row.names = FALSE, col.names = TRUE)

############################################################################################
############################################################################################
############################################################################################ SAMPLE6

SAMPLE6.file <- read.delim("SOX11_mutant2_anno_rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

SAMPLE6.df   <- data.frame(matrix(ncol = 4, nrow = dim(SAMPLE6.file)[1]))

colnames(SAMPLE6.df) <- c(paste(SAMPLE6, "gene", sep="."), 
                          paste(SAMPLE6, "count", sep="."), 
                          paste(SAMPLE6, "TPM", sep="."), 
                          paste(SAMPLE6, "FPKM", sep="."))

SAMPLE6.df[,1] = SAMPLE6.file$gene_id
SAMPLE6.df[,2] = SAMPLE6.file$expected_count
SAMPLE6.df[,3] = SAMPLE6.file$TPM
SAMPLE6.df[,4] = SAMPLE6.file$FPKM
                                 
dim(SAMPLE6.df)

# head(SAMPLE6.df)
# write.table(SAMPLE6.df, file=paste(SAMPLE6, "output.txt", sep="."), 
#                         quote = FALSE, sep = "\t", 
#                         row.names = FALSE, col.names = TRUE)

############################################################################################
############################################################################################ TO INTEGRATE THESE SAMPLES

### to integrate these dataframes : 
### SAMPLE1.df
### SAMPLE2.df
### SAMPLE3.df
### SAMPLE4.df
### SAMPLE5.df
### SAMPLE6.df

SAMPLE1.dt  <-  as.data.table(SAMPLE1.df)                            
SAMPLE2.dt  <-  as.data.table(SAMPLE2.df)                            
SAMPLE3.dt  <-  as.data.table(SAMPLE3.df)                            
SAMPLE4.dt  <-  as.data.table(SAMPLE4.df)                            
SAMPLE5.dt  <-  as.data.table(SAMPLE5.df)
SAMPLE6.dt  <-  as.data.table(SAMPLE6.df)

############################################################################################
############################################################################################
############################################################################################
############################################################################################

### in order to integrate the files with DATA.TABLE package : 

objects()
library(data.table)

setkeyv(genes.dt, c('gene_id'))

setkeyv(SAMPLE1.dt, c(paste(SAMPLE1, "gene", sep=".")) )
setkeyv(SAMPLE2.dt, c(paste(SAMPLE2, "gene", sep=".")) )
setkeyv(SAMPLE3.dt, c(paste(SAMPLE3, "gene", sep=".")) )
setkeyv(SAMPLE4.dt, c(paste(SAMPLE4, "gene", sep=".")) )
setkeyv(SAMPLE5.dt, c(paste(SAMPLE5, "gene", sep=".")) )
setkeyv(SAMPLE6.dt, c(paste(SAMPLE6, "gene", sep=".")) )

###############################################################################################
###############################################################################################

expression.all.samples <- genes.dt[SAMPLE1.dt,][SAMPLE2.dt,][SAMPLE3.dt,][SAMPLE4.dt,][SAMPLE5.dt,][SAMPLE6.dt,]

head(expression.all.samples)
dim(expression.all.samples)

# to write the integrated data file :

write.table(expression.all.samples, file=paste(name, ".INTEGRATED.file.ALL.samples.txt", sep=""),
                                    sep="\t", quote=FALSE,
                                    row.names = FALSE, col.names = TRUE)

# colnames(expression.all.samples)
# [1] "gene_id"           "gene_name"         "gene_type"        
# [4] "gene_status"       "transcript_id"     "transcript_type"  
# [7] "transcript_status" "transcript_name"   "Control_1.count"  
#[10] "Control_1.TPM"     "Control_1.FPKM"    "Control_2.count"  
#[13] "Control_2.TPM"     "Control_2.FPKM"    "Sox11_1.count"    
#[16] "Sox11_1.TPM"       "Sox11_1.FPKM"      "Sox11_2.count"    
#[19] "Sox11_2.TPM"       "Sox11_2.FPKM"      "Sox11MUT_1.count" 
#[22] "Sox11MUT_1.TPM"    "Sox11MUT_1.FPKM"   "Sox11MUT_2.count" 
#[25] "Sox11MUT_2.TPM"    "Sox11MUT_2.FPKM"

###############################################################################################
###############################################################################################
############################################################################################### no NA

# to eliminate the fields that are NA 

expression.all.samples.no.na.genes <- expression.all.samples[!is.na(expression.all.samples$gene_name),]

dim(expression.all.samples.no.na.genes)  # [1] 43280    26

###############################################################################################
############################################################################################### no MULTIPLICATES
###############################################################################################

# in order to follow the R code for DE, we would need to have UNIQUE GENE NAMES

length(expression.all.samples.no.na.genes$gene_name)
length(unique(expression.all.samples.no.na.genes$gene_name))

expression.all.samples.unique.genes <- expression.all.samples.no.na.genes[!duplicated(expression.all.samples.no.na.genes[,"gene_name"]),]

dim(expression.all.samples.unique.genes) # [1] 43280    26

write.table(expression.all.samples.unique.genes, 
                                    file=paste(name, ".INTEGRATED.file.ALL.samples.unique.gene.names.txt", sep=""),
                                    sep="\t", quote=FALSE,
                                    row.names = FALSE, col.names = TRUE)

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

rownames(expression.all.samples.unique.genes) <- expression.all.samples.unique.genes$gene_name

colnames(expression.all.samples.unique.genes)
# [1] "gene_id"           "gene_name"         "gene_type"        
# [4] "gene_status"       "transcript_id"     "transcript_type"  
# [7] "transcript_status" "transcript_name"   "Control_1.count"  
#[10] "Control_1.TPM"     "Control_1.FPKM"    "Control_2.count"  
#[13] "Control_2.TPM"     "Control_2.FPKM"    "Sox11_1.count"    
#[16] "Sox11_1.TPM"       "Sox11_1.FPKM"      "Sox11_2.count"    
#[19] "Sox11_2.TPM"       "Sox11_2.FPKM"      "Sox11MUT_1.count" 
#[22] "Sox11MUT_1.TPM"    "Sox11MUT_1.FPKM"   "Sox11MUT_2.count" 
#[25] "Sox11MUT_2.TPM"    "Sox11MUT_2.FPKM"  

################################################################################################################################
################################################################################################################################
############################################ shall it be useful to have here another dataframe : .. DF
################################################################################################################################

expression.all.samples.unique.genes.df <- as.data.frame(expression.all.samples.unique.genes)

rownames(expression.all.samples.unique.genes.df) <- expression.all.samples.unique.genes.df$gene_name

colnames(expression.all.samples.unique.genes.df)
# [1] "gene_id"           "gene_name"         "gene_type"        
# [4] "gene_status"       "transcript_id"     "transcript_type"  
# [7] "transcript_status" "transcript_name"   "Control_1.count"  
#[10] "Control_1.TPM"     "Control_1.FPKM"    "Control_2.count"  
#[13] "Control_2.TPM"     "Control_2.FPKM"    "Sox11_1.count"    
#[16] "Sox11_1.TPM"       "Sox11_1.FPKM"      "Sox11_2.count"    
#[19] "Sox11_2.TPM"       "Sox11_2.FPKM"      "Sox11MUT_1.count" 
#[22] "Sox11MUT_1.TPM"    "Sox11MUT_1.FPKM"   "Sox11MUT_2.count" 
#[25] "Sox11MUT_2.TPM"    "Sox11MUT_2.FPKM" 

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

### the following ID and GROUPS of SAMPLES :

# "gene_name"                 
# "Control_1.count"  
# "Control_2.count"  
# "Sox11_1.count"    
# "Sox11_2.count"    
# "Sox11MUT_1.count" 
# "Sox11MUT_2.count" 

### the NAMES of the COLUMNS are the following :

SAMPLE1.count <- paste(SAMPLE1, "count", sep=".")
SAMPLE2.count <- paste(SAMPLE2, "count", sep=".")
SAMPLE3.count <- paste(SAMPLE3, "count", sep=".")
SAMPLE4.count <- paste(SAMPLE4, "count", sep=".")
SAMPLE5.count <- paste(SAMPLE5, "count", sep=".")
SAMPLE6.count <- paste(SAMPLE6, "count", sep=".")

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ A NOTE :
### in order to access SPECIFIC COLUMNS :
### https://stackoverflow.com/questions/12603890/pass-column-name-in-data-table-using-variable
### https://stackoverflow.com/questions/18222286/dynamically-select-data-frame-columns-using-and-a-vector-of-column-names

# expression.all.samples.unique.genes[, get(SAMPLE1.count)]
# expression.all.samples.unique.genes[[(SAMPLE1.count)]]
# expression.all.samples.unique.genes.df[[(SAMPLE1.count)]]

# expression.all.samples.unique.genes[, get(SAMPLE2.count)]
# expression.all.samples.unique.genes[[(SAMPLE2.count)]]
# expression.all.samples.unique.genes.df[[(SAMPLE2.count)]]

# expression.all.samples.unique.genes[, get(SAMPLE3.count)]
# expression.all.samples.unique.genes[[(SAMPLE3.count)]]
# expression.all.samples.unique.genes.df[[(SAMPLE3.count)]]

# expression.all.samples.unique.genes[, get(SAMPLE4.count)]
# expression.all.samples.unique.genes[[(SAMPLE4.count)]]
# expression.all.samples.unique.genes.df[[(SAMPLE4.count)]]

# expression.all.samples.unique.genes[, get(SAMPLE5.count)]
# expression.all.samples.unique.genes[[(SAMPLE5.count)]]
# expression.all.samples.unique.genes.df[[(SAMPLE5.count)]]

# expression.all.samples.unique.genes[, get(SAMPLE6.count)]
# expression.all.samples.unique.genes[[(SAMPLE6.count)]]
# expression.all.samples.unique.genes.df[[(SAMPLE6.count)]]

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ A NOTE 
### in the continuation of the script, we can use the actual COLUMN NAMES :

# "gene_name" 
# "Control_1.count", 
# "Control_2.count",  
# "Sox11_1.count", 
# "Sox11_2.count"                
# "Sox11MUT_1.count" 
# "Sox11MUT_2.count" 

################################################################################################################################ CORRELATION
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ A NOTE
### the CORRELATION COEFFICIENT between the REPLICATES :

cor(expression.all.samples.unique.genes.df$Control_1.count, 
    expression.all.samples.unique.genes.df$Control_2.count)   # 0.96


cor(expression.all.samples.unique.genes.df$Sox11_1.count,
    expression.all.samples.unique.genes.df$Sox11_2.count)     # 0.96


cor(expression.all.samples.unique.genes.df$Sox11MUT_1.count,
    expression.all.samples.unique.genes.df$Sox11MUT_2.count)  # 0.97


### although the correlation coefficients on the other combinations of samples also looks very high ..

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
################################################################################################################################ COMPARISON1

genes.counts.Control.Sox11 <- subset(expression.all.samples.unique.genes.df, select=c("Control_1.count", 
                                                                                      "Control_2.count",  
                                                                                      "Sox11_1.count", 
                                                                                      "Sox11_2.count"))

head(genes.counts.Control.Sox11)
dim(genes.counts.Control.Sox11)  

################################################################################################################################ 
################################################################################################################################
################################################################################################################################

eset <- genes.counts.Control.Sox11
eset_name <- deparse(substitute(genes.counts.Control.Sox11)) ### in order to get the name of the DF

################################################################################################################################
################################################################# GROUP and SUBJECT (ie REPLICATE)
################################################################################################################################

group <- factor(c("Control", "Control", "Sox11", "Sox11"))
subject <- factor(c(1,2,1,2))

################################################################################################################################ 
################################################################################################################################
################################################################################################################################
#################################################################
################################################################# compute the CPM and display the DENSITIES :
################################################################################################################################

cpm_before_filtering  <- cpm(eset)
lcpm_before_filtering <- cpm(eset, log = TRUE)

# To visualize the distribution of gene expression levels BEFORE FILTERING

png(paste(name, eset_name, "plot.DENSITIES.before.FILTERING.png", sep="."))
plotDensities(lcpm_before_filtering, legend = TRUE, main = "before filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################# making a DGEList :
################################################################################################################################

x <- DGEList(counts=eset, group=group)

################################################################################################################################
################################################################################################################################
################################################################# FILTERING :
################################################################################################################################
### filtering the genes based on CPM :

keep.exprs <- rowSums( cpm(x) > 1) >= 2
x <- x[keep.exprs, ,keep.lib.sizes = FALSE]

dim(x)
x$samples$lib.size <- colSums(x$counts)

################################################################################################################################
################################################################################################################################
################################################################# and displaying the DENSITIES :
################################################################################################################################

lcpm_after_filtering <- cpm(x, log=TRUE)

# Visualize the distribution of gene expression levels AFTER FILTERING

png(paste(name, eset_name, "plot.DENSITIES.after.FILTERING.png", sep="."))
plotDensities(lcpm_after_filtering, legend = FALSE, main = "after filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# DESIGN and CONTRAST MATRIX :
################################################################################################################################

design <- model.matrix(~0+group+subject)

colnames(design) <- gsub("group", "", colnames(design))

# contrast.matrix <- makeContrasts(ControlvsSox11 = Control-Sox11, levels=design)    #### here using the SPECIFIC NAMES of SAMPLES

  contrast.matrix <- makeContrasts(ControlvsSox11 = Sox11-Control, levels=design)    #### changing the ORDER of SAMPLES

# design
#  Control Sox11 subject2
#1       1     0        0
#2       1     0        1
#3       0     1        0
#        0     1        1

# contrast.matrix 
#          Contrasts
# Levels     ControlvsSox11
#  Control              -1
#  Sox11                 1
#  subject2              0

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# NORMALIZATION :
################################################################################################################################
# Let's calculate the normalization factors for our data

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# VOOM :
################################################################################################################################
### using the VOOM transformation :

v <- voom(x, design, plot=FALSE)

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# LINEAR FIT :
################################################################################################################################
### the LINEAR FIT in LIMMA :

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)

efit <- eBayes(vfit)

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# MEAN VARIANCE TREND :
################################################################################################################################
### displaying the MEAN-VARIANCE TREND : 

png(paste(name, eset_name, "plot.MEAN.VARIANCE.after.NORM.png", sep="."))
plotSA(efit, main = "Final model: Mean−variance trend")
dev.off()

# Tabulate the results

summary(decideTests(efit))

#       ControlvsSox11
#Down             1837
#NotSig          10193
#Up               1487

# If the magnitude of the effect size is important for your downstream analysis, 
# you can specify a minimal log-fold-change with the function "treat" instead of using "eBayes"

# tfit <- treat(vfit, lfc = 1)
# dt <- decideTests(tfit)
# summary(dt)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ PRINTING THE RESULTS : 
################################################################################################################################
### obtaining and writing the results :

results <- topTable(efit, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns
results$Gene <- rownames(results)

write.table(results, 
            file=paste(name, eset_name, "RESULTS.limma.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ SEPARATING the DEG : 
################################################################################################################################

### computing the DEG for FDR < 0.05 :

results.deg <- results[results$adj.P.Val < 0.05,]

### computing the DEG for FDR < 0.05 and FC > 1.2 : UP REGULATED GENES

results.deg.up  <- results[(results$adj.P.Val < 0.05) & 
                           (results$logFC > log2(1.2) ), ]

dim(results.deg.up) # 1487

write.table(results.deg.up, 
            file=paste(name, eset_name, "RESULTS.limma.UP.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)


### computing the DEG for FDR < 0.05 and FC > 1.2 : DOWN REGULATED GENES

results.deg.down <- results[(results$adj.P.Val < 0.05) & 
                            (results$logFC < -log2(1.2) ), ]

dim(results.deg.down) # 1837

write.table(results.deg.down, 
            file=paste(name, eset_name, "RESULTS.limma.DOWN.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ TO RECORD THE GENES for integration at a later :
################################################################################################################################

genes.counts.Control.Sox11.deg      <- results.deg 
genes.counts.Control.Sox11.deg.up   <- results.deg.up 
genes.counts.Control.Sox11.deg.down <- results.deg.down 

genes.counts.Control.Sox11.limma <- results

dim(genes.counts.Control.Sox11.deg)      # [1] 3324    
dim(genes.counts.Control.Sox11.deg.up)   # [1] 1487   
dim(genes.counts.Control.Sox11.deg.down) # [1] 1837    
dim(genes.counts.Control.Sox11.limma)    # [1] 13517

################################################################################################################################
################################################################################################################################
################################################################
################################################################ TO RENAME the COLUMNS :            
################################################################################################################################

colnames(genes.counts.Control.Sox11.limma)[colnames(genes.counts.Control.Sox11.limma)=="logFC"]     <- paste("logFC", "Control", "Sox11", sep=".")
colnames(genes.counts.Control.Sox11.limma)[colnames(genes.counts.Control.Sox11.limma)=="AveExpr"]   <- paste("AveExpr", "Control", "Sox11", sep=".")
colnames(genes.counts.Control.Sox11.limma)[colnames(genes.counts.Control.Sox11.limma)=="t"]         <- paste("t", "Control", "Sox11", sep=".")
colnames(genes.counts.Control.Sox11.limma)[colnames(genes.counts.Control.Sox11.limma)=="P.Value"]   <- paste("P.Value", "Control", "Sox11", sep=".")
colnames(genes.counts.Control.Sox11.limma)[colnames(genes.counts.Control.Sox11.limma)=="adj.P.Val"] <- paste("adj.P.Val", "Control", "Sox11", sep=".")
colnames(genes.counts.Control.Sox11.limma)[colnames(genes.counts.Control.Sox11.limma)=="B"]         <- paste("B", "Control", "Sox11", sep=".")

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
################################################################################################################################ COMPARISON2

genes.counts.Control.Sox11MUT <- subset(expression.all.samples.unique.genes.df, select=c("Control_1.count", 
                                                                                         "Control_2.count",  
                                                                                         "Sox11MUT_1.count", 
                                                                                         "Sox11MUT_2.count"))

head(genes.counts.Control.Sox11MUT)
dim(genes.counts.Control.Sox11MUT)

################################################################################################################################ LIMMA-VOOM
################################################################################################################################

eset <- genes.counts.Control.Sox11MUT
eset_name <- deparse(substitute(genes.counts.Control.Sox11MUT)) ### in order to get the name of the DF

################################################################################################################################
################################################################# GROUP and SUBJECT (ie REPLICATE)
################################################################################################################################

group <- factor(c("Control", "Control", "Sox11MUT", "Sox11MUT"))
subject <- factor(c(1,2,1,2))

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# compute the CPM and display the DENSITIES :
################################################################################################################################

cpm_before_filtering  <- cpm(eset)
lcpm_before_filtering <- cpm(eset, log = TRUE)

# Visualize the distribution of gene expression levels BEFORE FILTERING

png(paste(name, eset_name, "plot.DENSITIES.before.FILTERING.png", sep="."))
plotDensities(lcpm_before_filtering, legend = TRUE, main = "before filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################# making a DGEList :
################################################################################################################################

x <- DGEList(counts=eset, group=group)

################################################################################################################################
################################################################################################################################
################################################################# FILTERING :
################################################################################################################################
### filtering the genes based on CPM :

keep.exprs <- rowSums( cpm(x) > 1) >= 2
x <- x[keep.exprs, ,keep.lib.sizes = FALSE]

dim(x)
x$samples$lib.size <- colSums(x$counts)

################################################################################################################################
################################################################################################################################
################################################################# and displaying the DENSITIES :
################################################################################################################################

lcpm_after_filtering <- cpm(x, log=TRUE)

# Visualize the distribution of gene expression levels AFTER FILTERING

png(paste(name, eset_name, "plot.DENSITIES.after.FILTERING.png", sep="."))
plotDensities(lcpm_after_filtering, legend = FALSE, main = "after filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# DESIGN and CONTRAST MATRIX :
################################################################################################################################

design <- model.matrix(~0+group+subject)

colnames(design) <- gsub("group", "", colnames(design))

# contrast.matrix <- makeContrasts(ControlvsSox11MUT = Control-Sox11MUT, levels=design)    #### here using the SPECIFIC NAMES of SAMPLES

  contrast.matrix <- makeContrasts(ControlvsSox11MUT = Sox11MUT-Control, levels=design)    #### changing the ORDER of the SAMPLES

design
#  Control Sox11MUT subject2
#1       1        0        0
#2       1        0        1
#3       0        1        0
#4       0        1        1

contrast.matrix
#          Contrasts
# Levels     ControlvsSox11
#  Control              -1
#  Sox11MUT              1
# subject2               0

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# NORMALIZATION :
################################################################################################################################
# Let's calculate the normalization factors for our data

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# VOOM :
################################################################################################################################
### using the VOOM transformation :

v <- voom(x, design, plot=FALSE)

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# LINEAR FIT :
################################################################################################################################
### the LINEAR FIT in LIMMA :

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)

efit <- eBayes(vfit)

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# MEAN VARIANCE TREND :
################################################################################################################################
### displaying the MEAN-VARIANCE TREND : 

png(paste(name, eset_name, "plot.MEAN.VARIANCE.after.NORM.png", sep="."))
plotSA(efit, main = "Final model: Mean−variance trend")
dev.off()

# Tabulate the results

summary(decideTests(efit))
#       ControlvsSox11
#Down              902
#NotSig          11725
#Up                911

# If the magnitude of the effect size is important for your downstream analysis, 
# you can specify a minimal log-fold-change with the function "treat" instead of using "eBayes"

# tfit <- treat(vfit, lfc = 1)
# dt <- decideTests(tfit)
# summary(dt)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ PRINTING THE RESULTS : 
################################################################################################################################

### obtaining and writing the results :
results <- topTable(efit, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns
results$Gene <- rownames(results)

write.table(results, 
            file=paste(name, eset_name, "RESULTS.limma.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ SEPARATING the DEG : 
################################################################################################################################

### computing the DEG for FDR < 0.05 :

results.deg <- results[results$adj.P.Val < 0.05,]

### computing the DEG for FDR < 0.05 and FC > 1.2 : UP REGULATED GENES

results.deg.up  <- results[(results$adj.P.Val < 0.05) & 
                           (results$logFC > log2(1.2) ), ]

dim(results.deg.up) # 911

write.table(results.deg.up, 
            file=paste(name, eset_name, "RESULTS.limma.UP.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)


### computing the DEG for FDR < 0.05 and FC > 1.2 : DOWN REGULATED GENES

results.deg.down <- results[(results$adj.P.Val < 0.05) & 
                            (results$logFC < -log2(1.2) ), ]

dim(results.deg.down) # 902

write.table(results.deg.down, 
            file=paste(name, eset_name, "RESULTS.limma.DOWN.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ TO RECORD THE GENES :
################################################################################################################################

genes.counts.Control.Sox11MUT.deg      <- results.deg 
genes.counts.Control.Sox11MUT.deg.up   <- results.deg.up 
genes.counts.Control.Sox11MUT.deg.down <- results.deg.down 

genes.counts.Control.Sox11MUT.limma <- results

dim(genes.counts.Control.Sox11MUT.deg)      # 1813
dim(genes.counts.Control.Sox11MUT.deg.up)   # 911
dim(genes.counts.Control.Sox11MUT.deg.down) # 902  
dim(genes.counts.Control.Sox11MUT.limma)    # 13538

################################################################################################################################
################################################################################################################################
################################################################
################################################################ TO RENAME the COLUMNS :            
################################################################################################################################

colnames(genes.counts.Control.Sox11MUT.limma)[colnames(genes.counts.Control.Sox11MUT.limma)=="logFC"]     <- paste("logFC", "Control", "Sox11MUT", sep=".")
colnames(genes.counts.Control.Sox11MUT.limma)[colnames(genes.counts.Control.Sox11MUT.limma)=="AveExpr"]   <- paste("AveExpr", "Control", "Sox11MUT", sep=".")
colnames(genes.counts.Control.Sox11MUT.limma)[colnames(genes.counts.Control.Sox11MUT.limma)=="t"]         <- paste("t", "Control", "Sox11MUT", sep=".")
colnames(genes.counts.Control.Sox11MUT.limma)[colnames(genes.counts.Control.Sox11MUT.limma)=="P.Value"]   <- paste("P.Value", "Control", "Sox11MUT", sep=".")
colnames(genes.counts.Control.Sox11MUT.limma)[colnames(genes.counts.Control.Sox11MUT.limma)=="adj.P.Val"] <- paste("adj.P.Val", "Control", "Sox11MUT", sep=".")
colnames(genes.counts.Control.Sox11MUT.limma)[colnames(genes.counts.Control.Sox11MUT.limma)=="B"]         <- paste("B", "Control", "Sox11MUT", sep=".")

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
################################################################################################################################ INTEGRATE all these 3 dataframes :

# expression.all.samples.unique.genes (or : expression.all.samples.unique.genes.df ?)
# genes.counts.Control.Sox11.limma 
# genes.counts.Control.Sox11MUT.limma 

# genes.counts.Control.Sox11.limma.dt     <-  as.data.table(genes.counts.Control.Sox11.limma)
# genes.counts.Control.Sox11MUT.limma.dt  <-  as.data.table(genes.counts.Control.Sox11MUT.limma)


################################################################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################ in order to integrate the files with DATA.TABLE : 
################################################################################################################################

# objects()
# library(data.table)

# setkeyv(expression.all.samples.unique.genes, c('gene_name'))

# setkeyv(genes.counts.Control.Sox11.limma.dt,    c('Gene'))
# setkeyv(genes.counts.Control.Sox11MUT.limma.dt, c('Gene')) 

################################################################################################################################
###############################################################################################
###############################################################################################
################################################################################################################################

# expression.all.samples.and.DEG <- expression.all.samples.unique.genes[genes.counts.Control.Sox11.limma.dt,][genes.counts.Control.Sox11MUT.limma.dt,]

# head(expression.all.samples.and.DEG)
# dim(expression.all.samples.and.DEG)

# writing the integrated data file :
# write.table(expression.all.samples.and.DEG, file=paste(name, ".INTEGRATED.file.ALL.samples.and.SETS.DEG.txt", sep=""),
#                                            sep="\t", quote=FALSE,
#                                            row.names = FALSE, col.names = TRUE)

################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ INTEGRATE all these 3 dataframes :

# expression.all.samples.unique.genes (or : expression.all.samples.unique.genes.df ?)
# genes.counts.Control.Sox11.limma 
# genes.counts.Control.Sox11MUT.limma 


integration.step1 <- merge(expression.all.samples.unique.genes,
                           genes.counts.Control.Sox11.limma,
                           by.x = "gene_name" ,
                           by.y = "Gene",
                           all.x = TRUE)

dim(integration.step1) #  43280    27

################################################################################################################################
###############################################################################################
################################################################################################################################

integration.step2 <- merge(integration.step1,
                           genes.counts.Control.Sox11MUT.limma,
                           by.x = "gene_name" ,
                           by.y = "Gene",
                           all.x = TRUE)

dim(integration.step2) # 43280    31

# colnames(integration.step2) # 43280    31
# [1] "gene_name"                  "gene_id"                   
# [3] "gene_type"                  "gene_status"               
# [5] "transcript_id"              "transcript_type"           
# [7] "transcript_status"          "transcript_name"           
# [9] "Control_1.count"            "Control_1.TPM"             
#[11] "Control_1.FPKM"             "Control_2.count"           
#[13] "Control_2.TPM"              "Control_2.FPKM"            
#[15] "Sox11_1.count"              "Sox11_1.TPM"               
#[17] "Sox11_1.FPKM"               "Sox11_2.count"             
#[19] "Sox11_2.TPM"                "Sox11_2.FPKM"              
#[21] "Sox11MUT_1.count"           "Sox11MUT_1.TPM"            
#[23] "Sox11MUT_1.FPKM"            "Sox11MUT_2.count"          
#[25] "Sox11MUT_2.TPM"             "Sox11MUT_2.FPKM"           
#[27] "logFC.Control.Sox11"        "AveExpr.Control.Sox11"     
#[29] "t.Control.Sox11"            "P.Value.Control.Sox11"     
#[31] "adj.P.Val.Control.Sox11"    "B.Control.Sox11"           
#[33] "logFC.Control.Sox11MUT"     "AveExpr.Control.Sox11MUT"  
#[35] "t.Control.Sox11MUT"         "P.Value.Control.Sox11MUT"  
#[37] "adj.P.Val.Control.Sox11MUT" "B.Control.Sox11MUT"        

write.table(integration.step2, 
            file=paste(name, "INTEGRATED.file.ALL.samples.and.SETS.DEG.large.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

###############################################################################################################################
###############################################################################################################################

### for simplicity 
### z <- integration.step2

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
################################################################################################################################ re_compute DEG
################################################################################################################################ CTRL vs SOX11
################################################################################################################################ full info

### computing the DEG for FDR < 0.05 and FC > 1.2 : UP REGULATED GENES i.e. HIGHER in CONTROL

genes.deg.up.Sox11  <- integration.step2[(integration.step2$adj.P.Val.Control.Sox11 < 0.05) & 
                                         (integration.step2$logFC.Control.Sox11 > log2(1.2) ), ]

dim(genes.deg.up.Sox11) # 1487

write.table(genes.deg.up.Sox11, 
            file=paste(name, "z.RESULTS.limma.compare.CTRL.SOX11.genes.UP.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)


### computing the DEG for FDR < 0.05 and FC > 1.2 : DOWN REGULATED GENES i.e. LOWER in CONTROL

genes.deg.down.Sox11 <- integration.step2[(integration.step2$adj.P.Val.Control.Sox11  < 0.05) & 
                                          (integration.step2$logFC.Control.Sox11  < -log2(1.2) ), ]

dim(genes.deg.down.Sox11) # 1837

write.table(genes.deg.down.Sox11, 
            file=paste(name, "z.RESULTS.limma.compare.CTRL.SOX11.genes.DOWN.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

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
################################################################################################################################ re_compute DEG
################################################################################################################################ CTRL vs SOX11MUT
################################################################################################################################ full info

### computing the DEG for FDR < 0.05 and FC > 1.2 : UP REGULATED GENES i.e. HIGHER in CONTROL

genes.deg.up.Sox11MUT  <- integration.step2[(integration.step2$adj.P.Val.Control.Sox11MUT < 0.05) & 
                                            (integration.step2$logFC.Control.Sox11MUT > log2(1.2) ), ]

dim(genes.deg.up.Sox11MUT) # 911

write.table(genes.deg.up.Sox11MUT, 
            file=paste(name, "z.RESULTS.limma.compare.CTRL.SOX11MUT.genes.UP.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)


### computing the DEG for FDR < 0.05 and FC > 1.2 : DOWN REGULATED GENES i.e. LOWER in CONTROL

genes.deg.down.Sox11MUT <- integration.step2[(integration.step2$adj.P.Val.Control.Sox11MUT  < 0.05) & 
                                             (integration.step2$logFC.Control.Sox11MUT < -log2(1.2) ), ]

dim(genes.deg.down.Sox11MUT) # 902

write.table(genes.deg.down.Sox11MUT, 
            file=paste(name, "z.RESULTS.limma.compare.CTRL.SOX11MUT.genes.DOWN.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

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
################################################################################################################################ compare UP_genes

sets_up <- list( up_Sox11 = genes.deg.up.Sox11$gene_name,
                 up_Sox11MUT = genes.deg.up.Sox11MUT$gene_name) 

sets_up_venn <- venn.diagram(sets_up, filename=NULL)

png(paste(name, "z.RESULTS.limma.compare.genes.UP.venn.diagram.png", sep="."))
grid.newpage()
grid.draw(sets_up_venn)
dev.off()

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
################################################################################################################################ compare DOWN_GENES

sets_down <- list( down_Sox11 = genes.deg.down.Sox11$gene_name,
                   down_Sox11MUT = genes.deg.down.Sox11MUT$gene_name) 

sets_down_venn <- venn.diagram(sets_down, filename=NULL)

png(paste(name, "z.RESULTS.limma.compare.genes.DOWN.venn.diagram.png", sep="."))
grid.newpage()
grid.draw(sets_down_venn)
dev.off()

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
################################################################################################################################ compare both UP and DOWN

sets_all <- list( up_Sox11 = genes.deg.up.Sox11$gene_name,
                  up_Sox11MUT = genes.deg.up.Sox11MUT$gene_name, 
                  down_Sox11 = genes.deg.down.Sox11$gene_name,
                  down_Sox11MUT = genes.deg.down.Sox11MUT$gene_name) 

sets_all_venn <- venn.diagram(sets_all, filename=NULL)

png(paste(name, "z.RESULTS.limma.compare.genes.BOTH.venn.diagram.png", sep="."))
grid.newpage()
grid.draw(sets_all_venn)
dev.off()

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ code enrichR

### at the end, we are going to run a script that calls enrichR for each data set 

# library("enrichR")

# dbs <- c("GO_Biological_Process_2018",
#          "GO_Cellular_Component_2018",
#          "GO_Molecular_Function_2018",
#          "DSigDB",
#          "Genome_Browser_PWMs",
#          "TRANSFAC_and_JASPAR_PWMs",
#          "ENCODE_TF_ChIP-seq_2014",
#          "ENCODE_TF_ChIP-seq_2015",
#          "ChEA_2016",
#          "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
#          "KEGG_2016",
#          "WikiPathways_2016",
#          "Reactome_2016",
#          "BioCarta_2016",
#          "Panther_2016",
#          "NCI-Nature_2016",
#          "OMIM_Disease",
#          "OMIM_Expanded",
#          "MSigDB_Computational",
#          "MSigDB_Oncogenic_Signatures",
#          "Chromosome_Location")

# list_genes <- rownames(LIST_GENES)

# enriched <- enrichr(list_genes, dbs)

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
######################################################################################################
######################################################################################################
########## printing the SELECTED databases : the R  code from the previous DOCUMENT 
########## to use the script in order to perform GO and PATHWAYS enrichment on the sets of DEG
########## we will have to add the script that does GSEA on the gene datasets.
######################################################################################################
######################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ comparison3

### here to add a piece of R code that directly compares WT versus MUTANT :

# colnames(expression.all.samples.unique.genes.df)
# [1] "gene_id"           "gene_name"         "gene_type"        
# [4] "gene_status"       "transcript_id"     "transcript_type"  
# [7] "transcript_status" "transcript_name"   "Control_1.count"  
#[10] "Control_1.TPM"     "Control_1.FPKM"    "Control_2.count"  
#[13] "Control_2.TPM"     "Control_2.FPKM"    "Sox11_1.count"    
#[16] "Sox11_1.TPM"       "Sox11_1.FPKM"      "Sox11_2.count"    
#[19] "Sox11_2.TPM"       "Sox11_2.FPKM"      "Sox11MUT_1.count" 
#[22] "Sox11MUT_1.TPM"    "Sox11MUT_1.FPKM"   "Sox11MUT_2.count" 
#[25] "Sox11MUT_2.TPM"    "Sox11MUT_2.FPKM"  

### 've added a piece of R code below to the previous code, although the number of DE GENES between SOX11 and SOX11_MUT is 0 .. uff !! !!!!

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################ COMPARISON3

genes.counts.Sox11.Sox11MUT <- subset(expression.all.samples.unique.genes.df, select=c("Sox11_1.count", 
                                                                                       "Sox11_2.count",  
                                                                                       "Sox11MUT_1.count", 
                                                                                       "Sox11MUT_2.count"))

head(genes.counts.Sox11.Sox11MUT)
dim(genes.counts.Sox11.Sox11MUT)

################################################################################################################################ LIMMA-VOOM
################################################################################################################################

eset <- genes.counts.Sox11.Sox11MUT
eset_name <- deparse(substitute(genes.counts.Sox11.Sox11MUT)) ### in order to get the name of the DF

################################################################# GROUP and SUBJECT (ie REPLICATE)

group <- factor(c("Sox11", "Sox11", "Sox11MUT", "Sox11MUT"))
subject <- factor(c(1,2,1,2))

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# compute the CPM and display the DENSITIES :

cpm_before_filtering  <- cpm(eset)
lcpm_before_filtering <- cpm(eset, log = TRUE)

# Visualize the distribution of gene expression levels BEFORE FILTERING

png(paste(name, eset_name, "plot.DENSITIES.before.FILTERING.png", sep="."))
plotDensities(lcpm_before_filtering, legend = TRUE, main = "before filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
################################################################# making a DGEList :
################################################################################################################################

x <- DGEList(counts=eset, group=group)

################################################################################################################################
################################################################# FILTERING :
### filtering the genes based on CPM :
################################################################################################################################

keep.exprs <- rowSums( cpm(x) > 1) >= 2
x <- x[keep.exprs, ,keep.lib.sizes = FALSE]

dim(x)
x$samples$lib.size <- colSums(x$counts)

################################################################################################################################
################################################################################################################################
################################################################# and displaying the DENSITIES :
################################################################################################################################

lcpm_after_filtering <- cpm(x, log=TRUE)

# Visualize the distribution of gene expression levels AFTER FILTERING

png(paste(name, eset_name, "plot.DENSITIES.after.FILTERING.png", sep="."))
plotDensities(lcpm_after_filtering, legend = FALSE, main = "after filtering")
abline(v = 0, lty = 3)
dev.off()

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# DESIGN and CONTRAST MATRIX :
################################################################################################################################

design <- model.matrix(~0+group+subject)

colnames(design) <- gsub("group", "", colnames(design))

contrast.matrix <- makeContrasts(Sox11vsSox11MUT = Sox11MUT-Sox11, levels=design)    #### changing the ORDER of the SAMPLES

design
#  Sox11 Sox11MUT subject2
#1       1        0        0
#2       1        0        1
#3       0        1        0
#4       0        1        1

contrast.matrix
#             Contrasts
#  Levels     Sox11vsSox11MUT
#  Sox11                 -1
#  Sox11MUT               1
#  subject2               0

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# NORMALIZATION :
################################################################################################################################
# Let's calculate the normalization factors for our data

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# VOOM :
################################################################################################################################
### to use the VOOM transformation :

v <- voom(x, design, plot=FALSE)

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# LINEAR FIT :
################################################################################################################################
### the LINEAR FIT in LIMMA :

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)

efit <- eBayes(vfit)

################################################################################################################################ 
################################################################################################################################
#################################################################
################################################################# MEAN VARIANCE TREND :
################################################################################################################################
### to display the MEAN-VARIANCE TREND : 

png(paste(name, eset_name, "plot.MEAN.VARIANCE.after.NORM.png", sep="."))
plotSA(efit, main = "Final model: Mean−variance trend")
dev.off()

# Tabulate the results

summary(decideTests(efit))
#Down                 0
#NotSig           13355
#Up                   0

# If the magnitude of the effect size is important for your downstream analysis, 
# you can specify a minimal log-fold-change with the function "treat" instead of using "eBayes"

# tfit <- treat(vfit, lfc = 1)
# dt <- decideTests(tfit)
# summary(dt)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ PRINTING THE RESULTS : 
################################################################################################################################

### obtaining and writing the results :
results <- topTable(efit, coef=1, adjust="fdr", number=Inf)

### adding the rownames as columns
results$Gene <- rownames(results)

write.table(results, 
            file=paste(name, eset_name, "RESULTS.limma.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ SEPARATING the DEG : 
################################################################################################################################
### computing the DEG for FDR < 0.05 :

results.deg <- results[results$adj.P.Val < 0.05,]

### computing the DEG for FDR < 0.05 and FC > 1.2 : UP REGULATED GENES

results.deg.up  <- results[(results$adj.P.Val < 0.05) & 
                           (results$logFC > log2(1.2) ), ]

dim(results.deg.up) # 0

write.table(results.deg.up, 
            file=paste(name, eset_name, "RESULTS.limma.UP.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)


### computing the DEG for FDR < 0.05 and FC > 1.2 : DOWN REGULATED GENES

results.deg.down <- results[(results$adj.P.Val < 0.05) & 
                            (results$logFC < -log2(1.2) ), ]

dim(results.deg.down) # 0

write.table(results.deg.down, 
            file=paste(name, eset_name, "RESULTS.limma.DOWN.txt", sep="."),
            sep="\t", quote=FALSE,
            row.names = FALSE, col.names = TRUE)

################################################################################################################################ 
################################################################################################################################
################################################################
################################################################ TO RECORD THE GENES :
################################################################################################################################

genes.counts.Sox11.Sox11MUT.deg      <- results.deg 
genes.counts.Sox11.Sox11MUT.deg.up   <- results.deg.up 
genes.counts.Sox11.Sox11MUT.deg.down <- results.deg.down 

genes.counts.Sox11.Sox11MUT.limma <- results

dim(genes.counts.Sox11.Sox11MUT.deg)      # 
dim(genes.counts.Sox11.Sox11MUT.deg.up)   # 
dim(genes.counts.Sox11.Sox11MUT.deg.down) #   
dim(genes.counts.Sox11.Sox11MUT.limma)    # 

################################################################################################################################
################################################################################################################################
################################################################
################################################################ TO RENAME the COLUMNS :            
################################################################################################################################

colnames(genes.counts.Sox11.Sox11MUT.limma)[colnames(genes.counts.Sox11.Sox11MUT.limma)=="logFC"]     <- paste("logFC", "Sox11", "Sox11MUT", sep=".")
colnames(genes.counts.Sox11.Sox11MUT.limma)[colnames(genes.counts.Sox11.Sox11MUT.limma)=="AveExpr"]   <- paste("AveExpr", "Sox11", "Sox11MUT", sep=".")
colnames(genes.counts.Sox11.Sox11MUT.limma)[colnames(genes.counts.Sox11.Sox11MUT.limma)=="t"]         <- paste("t", "Sox11", "Sox11MUT", sep=".")
colnames(genes.counts.Sox11.Sox11MUT.limma)[colnames(genes.counts.Sox11.Sox11MUT.limma)=="P.Value"]   <- paste("P.Value", "Sox11", "Sox11MUT", sep=".")
colnames(genes.counts.Sox11.Sox11MUT.limma)[colnames(genes.counts.Sox11.Sox11MUT.limma)=="adj.P.Val"] <- paste("adj.P.Val", "Sox11", "Sox11MUT", sep=".")
colnames(genes.counts.Sox11.Sox11MUT.limma)[colnames(genes.counts.Sox11.Sox11MUT.limma)=="B"]         <- paste("B", "Sox11", "Sox11MUT", sep=".")

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
