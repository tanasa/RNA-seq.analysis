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

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
############################################################################ NAME="analysis.DEseq2"

ddssva = ddssvasva 
CORRECT = "SVA"

# ddssva <- DESeqDataSetFromMatrix( countData = round(as.matrix(counts.small)),
#                               colData=samples, 
#                               design = ~ group)

rowRanges(ddssva)
colData(ddssva)
assays(ddssva)
assay(ddssva)
length(rowRanges(ddssva))

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ DE ANALYSIS
################################################################################################################################
# computing the differential expression and extracting the CONTRASTS
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################

# ddssva <- DESeq(ddssva) 
# resultsNames(ddssva) 

# [1] "Intercept"              "group_lacSTZ_vs_lacCTR" "group_zfpCTR_vs_lacCTR"
# [4] "group_zfpSTZ_vs_lacCTR"

#> using pre-existing size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# res <- results(ddssva)

# head(res[order(res$pvalue),])
# tail(res[order(res$pvalue),])

#head(res[order(res$pvalue),])
#log2 fold change (MLE): group zfpSTZ vs lacCTR 
#Wald test p-value: group zfpSTZ vs lacCTR 
#DataFrame with 6 rows and 6 columns
#                     baseMean    log2FoldChange             lfcSE
#                    <numeric>         <numeric>         <numeric>
#Ddx3          1313.1303907189  9.02193651818142 0.836029219226851
#LOC100911498 1733.87019023218 -4.83960418567183 0.461054400769813
#Eif2s3y      2197.82733524275  14.5148361239418  1.49819308901979
#Fosl1        196.032857985014 -7.22399656906971  1.05013770635197
#Rgs1         169.474355049354  6.27315728776023   1.0047093895629
#Sectm1b      26.4652548079873  19.4095217311205  3.12976627475903
#                          stat               pvalue                 padj
#                     <numeric>            <numeric>            <numeric>
#Ddx3          10.7914129203819 3.77933841614821e-27 4.88252729982188e-23
#LOC100911498  -10.496818114286 8.93409710557449e-26 5.77098002534584e-22
#Eif2s3y       9.68822792624033 3.38348048440369e-22 1.45703947926704e-18
# tail(res[order(res$pvalue),])
#log2 fold change (MLE): group zfpSTZ vs lacCTR 
#Wald test p-value: group zfpSTZ vs lacCTR 
#DataFrame with 6 rows and 6 columns
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
################################################################################################################################ the CONTRASTS of INTEREST

lacCTR.vs.lacSTZ <- as.data.frame(results(ddssva, contrast=c("group", "lacSTZ", "lacCTR")))
head(lacCTR.vs.lacSTZ)
dim(lacCTR.vs.lacSTZ)
mcols(results(ddssva, contrast=c("group", "lacSTZ", "lacCTR")), use.names = TRUE)
lacCTR.vs.lacSTZ$GENE = rownames(lacCTR.vs.lacSTZ)

zfpCTR.vs.zfpSTZ <- as.data.frame(results(ddssva, contrast=c("group", "zfpSTZ", "zfpCTR")))
head(zfpCTR.vs.zfpSTZ)
dim(zfpCTR.vs.zfpSTZ)
mcols(results(ddssva, contrast=c("group", "zfpSTZ", "zfpCTR")), use.names = TRUE)
zfpCTR.vs.zfpSTZ$GENE = rownames(zfpCTR.vs.zfpSTZ)

lacCTR.vs.zfpCTR <- as.data.frame(results(ddssva, contrast=c("group", "zfpCTR", "lacCTR")))
head(lacCTR.vs.zfpCTR)
dim(lacCTR.vs.zfpCTR)
mcols(results(ddssva, contrast=c("group", "zfpCTR", "lacCTR")), use.names = TRUE)
lacCTR.vs.zfpCTR$GENE = rownames(lacCTR.vs.zfpCTR)

lacSTZ.vs.zfpSTZ <- as.data.frame(results(ddssva, contrast=c("group", "zfpSTZ", "lacSTZ")))
head(lacSTZ.vs.zfpSTZ)
dim(lacSTZ.vs.zfpSTZ)
mcols(results(ddssva, contrast=c("group", "zfpSTZ", "lacSTZ")), use.names = TRUE)
lacSTZ.vs.zfpSTZ$GENE = rownames(lacSTZ.vs.zfpSTZ)

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ to print the RESULTS in FILES 

write.table(as.data.frame(lacCTR.vs.lacSTZ),
            file=paste(NAME, "differential.expression.lacCTR.vs.lacSTZ", CORRECT, "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(zfpCTR.vs.zfpSTZ),
            file=paste(NAME, "differential.expression.zfpCTR.vs.zfpSTZ", CORRECT, "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacCTR.vs.zfpCTR),
            file=paste(NAME, "differential.expression.lacCTR.vs.zfpCTR", CORRECT, "results.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacSTZ.vs.zfpSTZ),
            file=paste(NAME, "differential.expression.lacSTZ.vs.zfpSTZ", CORRECT, "results.txt", sep="."), 
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

# dim(lacCTR.vs.lacSTZ)
# dim(zfpCTR.vs.zfpSTZ)
# dim(lacCTR.vs.zfpCTR)
# dim(lacSTZ.vs.zfpSTZ)

dim(subset(lacCTR.vs.lacSTZ, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
dim(subset(lacCTR.vs.lacSTZ, (padj < 0.1) &  (abs(log2FoldChange) > 2)))

dim(subset(zfpCTR.vs.zfpSTZ, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
dim(subset(zfpCTR.vs.zfpSTZ, (padj < 0.1) &  (abs(log2FoldChange) > 2)))
 
dim(subset(lacCTR.vs.zfpCTR, (pvalue < 0.05) &  (abs(log2FoldChange) > 2)))
dim(subset(lacCTR.vs.zfpCTR, (padj < 0.1) & (abs(log2FoldChange) > 2)))

dim(subset(lacSTZ.vs.zfpSTZ,  (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
dim(subset(lacSTZ.vs.zfpSTZ,  (padj < 0.1) & (abs(log2FoldChange) > 2)))

colnames(lacCTR.vs.lacSTZ)
colnames(zfpCTR.vs.zfpSTZ)
colnames(lacCTR.vs.zfpCTR)
colnames(lacSTZ.vs.zfpSTZ)

################################################################################################################################ 
################################################################################################################################

#> dim(subset(lacCTR.vs.lacSTZ, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
#[1] 425   7
#> dim(subset(lacCTR.vs.lacSTZ, (padj < 0.1) &  (abs(log2FoldChange) > 2)))
#[1] 96  7
#> 
#> dim(subset(zfpCTR.vs.zfpSTZ, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
#[1] 815   7
#> dim(subset(zfpCTR.vs.zfpSTZ, (padj < 0.1) &  (abs(log2FoldChange) > 2)))
#[1] 692   7
#>  
#> dim(subset(lacCTR.vs.zfpCTR, (pvalue < 0.05) &  (abs(log2FoldChange) > 2)))
#[1] 787   7
#> dim(subset(lacCTR.vs.zfpCTR, (padj < 0.1) & (abs(log2FoldChange) > 2)))
#[1] 573   7
#> 
#> dim(subset(lacSTZ.vs.zfpSTZ,  (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
#[1] 376   7
#> dim(subset(lacSTZ.vs.zfpSTZ,  (padj < 0.1) & (abs(log2FoldChange) > 2)))
#[1] 80  7

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 

resultsNames(ddssva)
resultsNames(ddssva)

lacCTR.vs.lacSTZ.lfc <- as.data.frame(lfcShrink(ddssva, contrast=c("group", "lacSTZ", "lacCTR")))
lacCTR.vs.lacSTZ.lfc$GENE = rownames(lacCTR.vs.lacSTZ.lfc)

zfpCTR.vs.zfpSTZ.lfc <- as.data.frame(lfcShrink(ddssva, contrast=c("group", "zfpSTZ", "zfpCTR")))
zfpCTR.vs.zfpSTZ.lfc$GENE = rownames(zfpCTR.vs.zfpSTZ.lfc)

lacCTR.vs.zfpCTR.lfc <- as.data.frame(lfcShrink(ddssva, contrast=c("group", "zfpCTR", "lacCTR")))
lacCTR.vs.zfpCTR.lfc$GENE = rownames(lacCTR.vs.zfpCTR.lfc)

lacSTZ.vs.zfpSTZ.lfc <- as.data.frame(lfcShrink(ddssva, contrast=c("group", "zfpSTZ", "lacSTZ")))
lacSTZ.vs.zfpSTZ.lfc$GENE = rownames(lacSTZ.vs.zfpSTZ.lfc)

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 

#write.table(as.data.frame(lacCTR.vs.lacSTZ.lfc),
#            file=paste(NAME, "differential.expression.lacCTR.vs.lacSTZ", CORRECT,  "results.LFC.txt", sep="."), 
#            quote = FALSE, 
#            row.names = FALSE,
#            col.names = TRUE, sep = "\t")

#write.table(as.data.frame(zfpCTR.vs.zfpSTZ.lfc),
#            file=paste(NAME, "differential.expression.zfpCTR.vs.zfpSTZ", CORRECT,  "results.LFC.txt", sep="."), 
#            quote = FALSE, 
#            row.names = FALSE,
#            col.names = TRUE, sep = "\t")

#write.table(as.data.frame(lacCTR.vs.zfpCTR.lfc),
#            file=paste(NAME, "differential.expression.lacCTR.vs.zfpCTR", CORRECT,  "results.LFC.txt", sep="."), 
#            quote = FALSE, 
#            row.names = FALSE,
#            col.names = TRUE, sep = "\t")

#write.table(as.data.frame(lacSTZ.vs.zfpSTZ.lfc),
#            file=paste(NAME, "differential.expression.lacSTZ.vs.zfpSTZ", CORRECT,  "results.LFC.txt", sep="."), 
#            quote = FALSE, 
#            row.names = FALSE,
#            col.names = TRUE, sep = "\t")

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 

# dim(lacCTR.vs.lacSTZ)
# dim(zfpCTR.vs.zfpSTZ)
# dim(lacCTR.vs.zfpCTR)
# dim(lacSTZ.vs.zfpSTZ)

dim(subset(lacCTR.vs.lacSTZ.lfc, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
dim(subset(lacCTR.vs.lacSTZ.lfc, (padj < 0.1) & (abs(log2FoldChange) > 2)))

dim(subset(zfpCTR.vs.zfpSTZ.lfc, (pvalue < 0.05)& (abs(log2FoldChange) > 2)))
dim(subset(zfpCTR.vs.zfpSTZ.lfc, (padj < 0.1) & (abs(log2FoldChange) > 2)))
 
dim(subset(lacCTR.vs.zfpCTR.lfc, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
dim(subset(lacCTR.vs.zfpCTR.lfc, (padj < 0.1) & (abs(log2FoldChange) > 2)))

dim(subset(lacSTZ.vs.zfpSTZ.lfc, (pvalue < 0.05) & (abs(log2FoldChange) > 2)))
dim(subset(lacSTZ.vs.zfpSTZ.lfc, (padj < 0.1) & (abs(log2FoldChange) > 2)))

################################################################################################################################ 
################################################################################################################################
################################################################################################################################ 
################################################################################################################################

write.table(as.data.frame(lacCTR.vs.lacSTZ.lfc),
            file=paste(NAME, "differential.expression.lacCTR.vs.lacSTZ", CORRECT, "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(zfpCTR.vs.zfpSTZ.lfc),
            file=paste(NAME, "differential.expression.zfpCTR.vs.zfpSTZ", CORRECT, "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacCTR.vs.zfpCTR.lfc),
            file=paste(NAME, "differential.expression.lacCTR.vs.zfpCTR", CORRECT, "results.lfc.txt", sep="."), 
            quote = FALSE, 
            row.names = FALSE,
            col.names = TRUE, sep = "\t")

write.table(as.data.frame(lacSTZ.vs.zfpSTZ.lfc),
            file=paste(NAME, "differential.expression.lacSTZ.vs.zfpSTZ", CORRECT, "results.lfc.txt", sep="."), 
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
################################################################################################################################
################################################################################################################################ 
################################################################################################################################
################################################################################################################################
################################################################################################################################ INTEGRATING these DATAFRAMES
################################################################################################################################
################################################################################################################################ INTEGRATING these FILES 

# using COUNTS
# counts = read.delim("KALLISTO.RefSeq.INTEGRATED.TPM.txt.no.dup.genes.txt", 
#                      sep="\t", header=T, stringsAsFactors=F)

# lacCTR.vs.lacSTZ <- as.data.frame(results(ddssva, contrast=c("group", "lacSTZ", "lacCTR")))
# zfpCTR.vs.zfpSTZ <- as.data.frame(results(ddssva, contrast=c("group", "zfpSTZ", "zfpCTR")))
# lacCTR.vs.zfpCTR <- as.data.frame(results(ddssva, contrast=c("group", "zfpCTR", "lacCTR")))
# lacSTZ.vs.zfpSTZ <- as.data.frame(results(ddssva, contrast=c("group", "zfpSTZ", "lacSTZ")))

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

png(paste(NAME, "part3.MA", "lacCTR.vs.lacSTZ", CORRECT, "png", sep="."))
plot(log2(lacCTR.vs.lacSTZ$baseMean), lacCTR.vs.lacSTZ$log2FoldChange)
lines(lowess(log2(lacCTR.vs.lacSTZ$baseMean), lacCTR.vs.lacSTZ$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

png(paste(NAME, "part3.MA", "zfpCTR.vs.zfpSTZ", CORRECT, "png", sep="."))
plot(log2(zfpCTR.vs.zfpSTZ$baseMean), zfpCTR.vs.zfpSTZ$log2FoldChange)
lines(lowess(log2(zfpCTR.vs.zfpSTZ$baseMean), zfpCTR.vs.zfpSTZ$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

png(paste(NAME, "part3.MA", "lacCTR.vs.zfpCTR", CORRECT, "png", sep="."))
plot(log2(lacCTR.vs.zfpCTR$baseMean), lacCTR.vs.zfpCTR$log2FoldChange)
lines(lowess(log2(lacCTR.vs.zfpCTR$baseMean), lacCTR.vs.zfpCTR$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

png(paste(NAME, "part3.MA", "lacSTZ.vs.zfpSTZ", CORRECT, "png", sep="."))
plot(log2(lacSTZ.vs.zfpSTZ$baseMean), lacSTZ.vs.zfpSTZ$log2FoldChange)
lines(lowess(log2(lacSTZ.vs.zfpSTZ$baseMean), lacSTZ.vs.zfpSTZ$log2FoldChange), col="blue") # lowess line (x,y)
dev.off()

colnames(counts)
colnames(lacCTR.vs.lacSTZ) 
colnames(zfpCTR.vs.zfpSTZ) 
colnames(lacCTR.vs.zfpCTR) 
colnames(lacSTZ.vs.zfpSTZ) 

#> colnames(counts)
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
#> colnames(lacCTR.vs.lacSTZ) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE"          
#> colnames(zfpCTR.vs.zfpSTZ) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE"          
#> colnames(lacCTR.vs.zfpCTR) 
#[1] "baseMean"       "log2FoldChange" "lfcSE"          "stat"          
#[5] "pvalue"         "padj"           "GENE"          
#> colnames(lacSTZ.vs.zfpSTZ) 
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

counts.SMALL.TPM = counts.TPM

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

##############################################################################################################################
##############################################################################################################################

head(lacCTR.vs.lacSTZ) 
head(zfpCTR.vs.zfpSTZ) 
head(lacCTR.vs.zfpCTR) 
head(lacSTZ.vs.zfpSTZ) 

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
            file=paste(NAME, "the.INTEGRATED.TABLES", CORRECT, "RESULTS.txt", sep="."), 
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
            file=paste(NAME, "the.INTEGRATED.TABLES", CORRECT, "RESULTS.LFC.txt", sep="."), 
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
#### COMBAT-seq might be another option
################################################################################################################################
################################################################################################################################
################################################################################################################################
