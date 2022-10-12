###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
## the script has been updated on OCT 7th, 2019
## and again on DEC 2nd, 2019
## we make the BARPLOTS of the CELL CYCLE GENES
###############################################################################################
###############################################################################################

library("ggplot2")
library("reshape2")
library("data.table")
library("tidyr")
library("plyr")
library("dplyr")
library("gplots")
library("pheatmap")
library("RColorBrewer")
library("scatterplot3d")
library("limma")
library("edgeR")
library("Glimma")
library("DESeq2")

################################################################################################
################################################################################################

NAME <- "the.analysis.results.03OCT2019"

#################################################################################################
######################### to upload the data where we did integrate all the files with ALL GENES
########################################################################### the results from RSEM
########################################################################## the results from LIMMA

genes.expression.large <- read.delim("analysis.LIMMA.integrating.all.samples.all.genes.in.4.STEPS.on.03oct2019.txt",
                     sep="\t", header=T, stringsAsFactors=F)

head(genes.expression.large) 
dim(genes.expression.large)

######################################################### we would have to make a special ROWNAME, 
################################################### as some genes are present in multiple isoforms 

genes.expression.large$ID <- rownames(genes.expression.large) 
genes.expression.large$GENE_NAME_ID <- paste(genes.expression.large$GENE_NAME, 
                                             genes.expression.large$ID, sep=":")

head(genes.expression.large)
dim(genes.expression.large)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
################################################## transforming the DATA FRAME into a DATA TABLE :

genes.expression.large.dt <- as.data.table(genes.expression.large)

head(genes.expression.large.dt) 
dim(genes.expression.large.dt)

############################################################################################
############################################################################################
############################################################################################
############################################################################################

###### selecting the following fields below in order to make
###### the PCA plots
###### the MDS plots 
###### the BOXPLOTS 
###### the SCATTER PLOTS
###### the VOLCANO PLOTS
###### the HEATMAPS

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
################################################################### DATAFRAME of GENES COUNTS :

genes.expression.large.counts <- subset(genes.expression.large, 
                                               select=c("GENE_NAME_ID", 
                                               "DMSO1_lane1.count", "DMSO1_lane2.count",
                                               "DMSO2_lane1.count", "DMSO2_lane2.count",
                                               "DMSO3_lane1.count", "DMSO3_lane2.count", 
                                               "Aph1.count", "Aph2.count", "Aph3.count", 
                                               "Aph_KH7_1.count","Aph_KH7_2.count","Aph_KH7_3.count",
                                               "KH7_1.count", "KH7_2.count", "KH7_3.count",
                                               "Noc_1.count", "Noc_2.count", "Noc_3.count"), 
                                                na.rm = TRUE)

rownames(genes.expression.large.counts) <- genes.expression.large.counts$GENE_NAME_ID 
genes.expression.large.counts <- genes.expression.large.counts[,-1]

head(genes.expression.large.counts)
dim(genes.expression.large.counts)

#############################################################################################
#############################################################################################
############################################################# making a DATAFRAME based on TPM :

genes.expression.large.tpm <- subset(genes.expression.large, 
                                     select=c("GENE_NAME_ID", 
                                     "DMSO1_lane1.TPM", "DMSO1_lane2.TPM",
                                     "DMSO2_lane1.TPM", "DMSO2_lane2.TPM",
                                     "DMSO3_lane1.TPM", "DMSO3_lane2.TPM", 
                                     "Aph1.TPM", "Aph2.TPM", "Aph3.TPM", 
                                     "Aph_KH7_1.TPM","Aph_KH7_2.TPM","Aph_KH7_3.TPM",
                                     "KH7_1.TPM", "KH7_2.TPM", "KH7_3.TPM",
                                     "Noc_1.TPM", "Noc_2.TPM", "Noc_3.TPM" ), 
                                      na.rm = TRUE)

rownames(genes.expression.large.tpm) <- genes.expression.large.tpm$GENE_NAME_ID 
genes.expression.large.tpm <- genes.expression.large.tpm[,-1]

head(genes.expression.large.tpm)
dim(genes.expression.large.tpm)

############################################################################################
########################################################## making a DATAFRAME based on FPKM :

genes.expression.large.fpkm <- subset(genes.expression.large, 
                                      select=c("GENE_NAME_ID", 
                                      "DMSO1_lane1.FPKM", "DMSO1_lane2.FPKM",
                                      "DMSO2_lane1.FPKM", "DMSO2_lane2.FPKM",
                                      "DMSO3_lane1.FPKM", "DMSO3_lane2.FPKM", 
                                      "Aph1.FPKM", "Aph2.FPKM", "Aph3.FPKM", 
                                      "Aph_KH7_1.FPKM","Aph_KH7_2.FPKM","Aph_KH7_3.FPKM",
                                      "KH7_1.FPKM", "KH7_2.FPKM", "KH7_3.FPKM",
                                      "Noc_1.FPKM", "Noc_2.FPKM", "Noc_3.FPKM" ), 
                                       na.rm = TRUE)

rownames(genes.expression.large.fpkm) <- genes.expression.large.fpkm$GENE_NAME_ID 
genes.expression.large.fpkm <- genes.expression.large.fpkm[,-1]

head(genes.expression.large.fpkm)
dim(genes.expression.large.fpkm)



################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
######################################### to upload the data where we did integrate all the files with ALL GENES
######################################### the results from RSEM
######################################### the results from LIMMA

genes.expression.small <- read.delim("analysis.LIMMA.integrating.all.samples.with.DATA.TABLE.on.03oct2019.txt",
                                      sep="\t", header=T, stringsAsFactors=F)

head(genes.expression.small) 
dim(genes.expression.small)

###################################################################### we would have to make a special ROWNAME, 
############################################################### as some genes are present in multiple isoforms. 

genes.expression.small$ID <- rownames(genes.expression.small) 
genes.expression.small$GENE_NAME_ID <- paste(genes.expression.small$GENE_NAME, 
                                             genes.expression.small$ID, sep=":")

head(genes.expression.small)
dim(genes.expression.small)

############################################################################################
############################################################################################
############################################################################################
############################################################################################

genes.expression.small.NA <- subset(genes.expression.small, is.na(CHR))
dim(genes.expression.small.NA)     ### 291  95

genes.expression.small.non.NA <- subset(genes.expression.small, !is.na(CHR))
dim(genes.expression.small.non.NA) ### 12956    95

genes.expression.small <- genes.expression.small.non.NA
dim(genes.expression.small)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################### transforming the DATA FRAME into a DATA TABLE :

genes.expression.small.dt <- as.data.table(genes.expression.small)

head(genes.expression.small.dt) 
dim(genes.expression.small.dt)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
########################################################################## making a DATAFRAME of GENES COUNTS :

genes.expression.small.counts <- subset(genes.expression.small, 
                                               select=c("GENE_NAME_ID", 
                                               "DMSO1_lane1.count", "DMSO1_lane2.count",
                                               "DMSO2_lane1.count", "DMSO2_lane2.count",
                                               "DMSO3_lane1.count", "DMSO3_lane2.count", 
                                               "Aph1.count", "Aph2.count", "Aph3.count", 
                                               "Aph_KH7_1.count","Aph_KH7_2.count","Aph_KH7_3.count",
                                               "KH7_1.count", "KH7_2.count", "KH7_3.count",
                                               "Noc_1.count", "Noc_2.count", "Noc_3.count"), 
                                                na.rm = TRUE)

rownames(genes.expression.small.counts) <- genes.expression.small.counts$GENE_NAME_ID 
genes.expression.small.counts <- genes.expression.small.counts[,-1]

head(genes.expression.small.counts)
dim(genes.expression.small.counts)

##################################################################################
##################################################################################
##################################################################################
##################################################################################
################################################ making a DATAFRAME based on TPM :

genes.expression.small.tpm <- subset(genes.expression.small, 
                        select=c("GENE_NAME_ID", 
                       "DMSO1_lane1.TPM", "DMSO1_lane2.TPM",
                       "DMSO2_lane1.TPM", "DMSO2_lane2.TPM",
                       "DMSO3_lane1.TPM", "DMSO3_lane2.TPM", 
                       "Aph1.TPM", "Aph2.TPM", "Aph3.TPM", 
                       "Aph_KH7_1.TPM","Aph_KH7_2.TPM","Aph_KH7_3.TPM",
                       "KH7_1.TPM", "KH7_2.TPM", "KH7_3.TPM",
                       "Noc_1.TPM", "Noc_2.TPM", "Noc_3.TPM" ), 
                        na.rm = TRUE)

rownames(genes.expression.small.tpm) <- genes.expression.small.tpm$GENE_NAME_ID 
genes.expression.small.tpm <- genes.expression.small.tpm[,-1]

head(genes.expression.small.tpm)
dim(genes.expression.small.tpm)

##################################################################################
##################################################################################
##################################################################### the MEDIAN : 

median(genes.expression.small.tpm[,1],na.rm=T)
median(genes.expression.small.tpm[,2],na.rm=T)
median(genes.expression.small.tpm[,3],na.rm=T)
median(genes.expression.small.tpm[,4],na.rm=T)
median(genes.expression.small.tpm[,5],na.rm=T)
median(genes.expression.small.tpm[,6],na.rm=T)
median(genes.expression.small.tpm[,7],na.rm=T)
median(genes.expression.small.tpm[,8],na.rm=T)
median(genes.expression.small.tpm[,9],na.rm=T)
median(genes.expression.small.tpm[,10],na.rm=T)
median(genes.expression.small.tpm[,11],na.rm=T)
median(genes.expression.small.tpm[,12],na.rm=T)
median(genes.expression.small.tpm[,13],na.rm=T)
median(genes.expression.small.tpm[,14],na.rm=T)
median(genes.expression.small.tpm[,15],na.rm=T)
median(genes.expression.small.tpm[,16],na.rm=T)
median(genes.expression.small.tpm[,17],na.rm=T)
median(genes.expression.small.tpm[,18],na.rm=T)

##################################################################################
####################################################### the BOXPLOTS for the genes 

pdf(paste(NAME, ".small.matrix.boxplot.TPM.pdf", sep="")) 
    par(las=2)
    par(mar=c(8,4,2,2))

    boxplot(genes.expression.small.tpm, 
            ylim=c(0,60), 
            col=c(rep("red",6), rep("orange",3), rep("green",3),
                            rep("blue",3), rep("violet",3)),
            ylab="TPM", 
            main="TPM values of ~13240 genes",
            cex.main=0.8, cex.lab=0.8)
dev.off()

##################################################################################
############################################# printing the FILE with TPM values :

write.table(genes.expression.small.tpm, 
            file=paste(NAME, ".small.file.TPM.txt", sep=""), 
            sep="\t") 

##################################################################################
##################################################################################
##################################################################################
##################################################################################
############################################### making a DATAFRAME based on FPKM :

genes.expression.small.fpkm <- subset(genes.expression.small, 
                        select=c("GENE_NAME_ID", 
                       "DMSO1_lane1.FPKM", "DMSO1_lane2.FPKM",
                       "DMSO2_lane1.FPKM", "DMSO2_lane2.FPKM",
                       "DMSO3_lane1.FPKM", "DMSO3_lane2.FPKM", 
                       "Aph1.FPKM", "Aph2.FPKM", "Aph3.FPKM", 
                       "Aph_KH7_1.FPKM","Aph_KH7_2.FPKM","Aph_KH7_3.FPKM",
                       "KH7_1.FPKM", "KH7_2.FPKM", "KH7_3.FPKM",
                       "Noc_1.FPKM", "Noc_2.FPKM", "Noc_3.FPKM" ), 
                        na.rm = TRUE)

rownames(genes.expression.small.fpkm) <- genes.expression.small.fpkm$GENE_NAME_ID 
genes.expression.small.fpkm <- genes.expression.small.fpkm[,-1]

head(genes.expression.small.fpkm)
dim(genes.expression.small.fpkm)

##################################################################################
##################################################################################
########################################################## to look at the MEDIAN : 

median(genes.expression.small.fpkm[,1],na.rm=T)
median(genes.expression.small.fpkm[,2],na.rm=T)
median(genes.expression.small.fpkm[,3],na.rm=T)
median(genes.expression.small.fpkm[,4],na.rm=T)
median(genes.expression.small.fpkm[,5],na.rm=T)
median(genes.expression.small.fpkm[,6],na.rm=T)
median(genes.expression.small.fpkm[,7],na.rm=T)
median(genes.expression.small.fpkm[,8],na.rm=T)
median(genes.expression.small.fpkm[,9],na.rm=T)
median(genes.expression.small.fpkm[,10],na.rm=T)
median(genes.expression.small.fpkm[,11],na.rm=T)
median(genes.expression.small.fpkm[,12],na.rm=T)
median(genes.expression.small.fpkm[,13],na.rm=T)
median(genes.expression.small.fpkm[,14],na.rm=T)
median(genes.expression.small.fpkm[,15],na.rm=T)
median(genes.expression.small.fpkm[,16],na.rm=T)
median(genes.expression.small.fpkm[,17],na.rm=T)
median(genes.expression.small.fpkm[,18],na.rm=T)

##################################################################################
############################################# making the BOXPLOTS for the genes :

pdf(paste(NAME, ".small.matrix.boxplot.FPKM.pdf", sep="")) 
    par(las=2)
    par(mar=c(8,4,2,2))

    boxplot(genes.expression.small.fpkm, 
            ylim=c(0,60), 
            col=c(rep("red",6), rep("orange",3), rep("green",3),
                            rep("blue",3), rep("violet",3)),
            ylab="FPKM", 
            main="FPKM values of ~13240 genes",
            cex.main=0.8, cex.lab=0.8)
dev.off()

##################################################################################
############################################ printing the FILE with FPKM values :

write.table(genes.expression.small.fpkm, 
            file=paste(NAME, ".small.file.FPKM.txt", sep=""), 
            sep="\t") 

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

###### selecting the following fields below in order to make
###### the PCA plots
###### the MDS plots 
###### the BOXPLOTS 
###### the HEATMAPS
###### the SCATTER PLOTS
###### the VOLCANO PLOTS

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

# In order to set up the color palette :
# col.rainbow <- rainbow(12)
# col.topo <- topo.colors(12)
# col.terrain <- terrain.colors(12)
# palette(col.rainbow)

######################################################################################
######################################################################################
###################################################################################### on the small matrix
###################################################################################### PCA ANALYSIS
###################################################################################### 

colnames(genes.expression.small.fpkm)

# [1]  "DMSO1_lane1.FPKM" "DMSO1_lane2.FPKM" "DMSO2_lane1.FPKM" "DMSO2_lane2.FPKM"
# [5]  "DMSO3_lane1.FPKM" "DMSO3_lane2.FPKM" 
#      "Aph1.FPKM"        "Aph2.FPKM"       
# [9]  "Aph3.FPKM"        
#      "Aph_KH7_1.FPKM"   "Aph_KH7_2.FPKM"   "Aph_KH7_3.FPKM"  
# [13] "KH7_1.FPKM"      "KH7_2.FPKM"       "KH7_3.FPKM"       
#      "Noc_1.FPKM"      
# [17] "Noc_2.FPKM"      "Noc_3.FPKM"

group <- factor( c(rep("DMSO",6), rep("Aph", 3), rep("Aph_KH7", 3), 
                   rep("KH7", 3), rep("Noc", 3)) )

######################################################################################
######################################################################################
######################################################################################
###################################################################################### 

library(scatterplot3d)

pca <- prcomp(t(genes.expression.small.fpkm))

######################################################################################
######################################################################################
################################################################### plotting the pca$x  

colors = c(rep("red",6), rep("orange",3), rep("green",3), rep("blue",3), rep("violet",3))
                     
######################################################################################
######################################################################################
################################################################# printing the image :

pdf(paste(NAME, ".small.file.PCA.display.in.3D.pdf", sep=""), width=10, height=10)
       
s3d <- scatterplot3d(pca$x[,1:3], 
                    color=colors, 
                    pch=18,
                    main="PCA analysis of DEG", 
                    grid=TRUE, 
                    box=TRUE)

legend(# s3d$xyz.convert(7, 4, 4), 
       "left",
       legend = row.names(pca$x[,1:3]),
       col = colors,                                     ### it has to use COL, and not COLORS
       pch = 16)

#legend("right", 
#       legend = row.names(pca$x[,1:3]),
#       col = colors,                                    ### it has to use COL, and not COLORS
#       pch = 16)

#legend("right", 
#       legend = row.names(pca$x[,1:3]),
#       col = colors,                                    ### it has to use COL, and not COLORS
#       pch = 16, 
#       inset = -0.15, xpd = TRUE, horiz = TRUE)

dev.off()

######################################################################################
######################################################################################
######################################################### to get some inspiration from :
# http://www.sthda.com/english/wiki/scatterplot3d-3d-graphics-r-software-and-data-visualization

# legend(s3D$xyz.convert(7.5, 3, 4.5), 
#       legend = row.names(pca$x[,1:3]),
#       color = c(rep("red",6), rep("orange",3), rep("green",3),
#                                     rep("blue",3), rep("violet",3)), 
#       pch = 16)

######################################################################################
######################################################################################
######################################################################################
######################################################################################

pca.df <- data.frame(PCA1=pca$x[,1], 
                     PCA2=pca$x[,2], 
                     PCA3=pca$x[,3], 
                     group=group)

######################################################################################
######################################################################################
############################### we are plotting PCA1 vs PCA2

pdf(paste(NAME, ".small.file.PCA.display.PC1.vs.PC2.pdf", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

######################################################################################
######################################################################################
################################ we are plotting PCA1 vs PCA3

pdf(paste(NAME, ".small.file.PCA.display.PC1.vs.PC3.pdf", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA3, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC3") + 
        ggtitle("PCA analysis : PC1 vs PC3")
dev.off()

######################################################################################
######################################################################################
################################# we are plotting PCA2 vs PCA3

pdf(paste(NAME, ".small.file.PCA.display.PC2.vs.PC3.pdf", sep=""))
ggplot(pca.df, 
       aes(x=PCA2, y=PCA3, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC2", y="PC3") + 
        ggtitle("PCA analysis : PC2 vs PC3")
dev.off()

########################################################################################################################
########################################################################################################################
######################################################################################################################## MDS ANALYSIS
######################################################################################################################## on small matrix
######################################################################################################################## 

group <- factor( c(rep("DMSO",6), rep("Aph", 3), rep("Aph_KH7", 3), 
                   rep("KH7", 3), rep("Noc", 3)) )

### We can use the function plotMDS from LIMMA or we can use the function cmdscale :
### mds <- plotMDS(genes.expression.small.fpkm)
### mds.df <- data.frame(MDSx=mds$x, MDSy=mds$y, group=group)

mds <- cmdscale(dist(t(genes.expression.small.fpkm)))
mds.df <- data.frame(MDSx=mds[,1], MDSy=mds[,2], group=group)

### plot(cmdscale(dist(t(genes.expression.small.fpkm))))
### text(cmdscale(dist(t(genes.expression.small.fpkm))), 
###               labels=colnames(genes.expression.small.fpkm))

pdf(paste(NAME, ".MDS.display.MDS1.vs.MDS2.pdf", sep=""))
ggplot(mds.df, aes(x=MDSx, 
                   y=MDSy, 
                   color=group, 
                   label=rownames(mds.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", legend.title=element_blank(), 
                                    legend.key = element_blank()) +
       labs(x="MDS dimension 1", y="MDS dimension 2") +
       ggtitle("MDS display")
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################### the BOXPLOTS :

###### before/after NORMALIZATION 
###### NOT-NORMALIZED counts
###### NORMALIZED counts

########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################### the HEATMAPS :

# par(las=2)
# par(mar=c(8,4,2,2))
# par(cex.main=0.6)

# heatmap.2(as.matrix(genes.expression.small.fpkm), col=bluered(149),
#                    scale="row",trace="none",
#                    cexRow=0.6, cexCol=0.6, cex.main=0.6,
#                    Rowv=FALSE, symkey=FALSE, labRow=NA,
#                    key=T, keysize=1.5, density.info="none",
#                    main="heatmap of ~ 1300 genes")

########################################################################################################################
########################################################################################################################
#################################################################################### to look at the following columns : 

colnames(genes.expression.small)
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
#[61] "Noc_3.FPKM"        "ID"                "GENE_NAME_ID"     
#[64] "logFC.Aph"         "AveExpr.Aph"       "t.Aph"            
#[67] "P.Value.Aph"       "adj.P.Val.Aph"     "B.Aph"            
#[70] "Gene.Aph"          "ID.Aph"            "logFC.Aph_KH7"    
#[73] "AveExpr.Aph_KH7"   "t.Aph_KH7"         "P.Value.Aph_KH7"  
#[76] "adj.P.Val.Aph_KH7" "B.Aph_KH7"         "Gene.Aph_KH7"     
#[79] "ID.Aph_KH7"        "logFC.KH7"         "AveExpr.KH7"      
#[82] "t.KH7"             "P.Value.KH7"       "adj.P.Val.KH7"    
#[85] "B.KH7"             "Gene.KH7"          "ID.KH7"           
#[88] "logFC.Noc"         "AveExpr.Noc"       "t.Noc"            
#[91] "P.Value.Noc"       "adj.P.Val.Noc"     "B.Noc"            
#[94] "Gene.Noc"          "ID.Noc"  

########################################################################################################################
################################################ to select the following columns : based on a THRESHOLD :

## "adj.P.Val.Aph"      
## "adj.P.Val.Aph_KH7" 
## "adj.P.Val.KH7"    
## "adj.P.Val.Noc"

FDR_THRESHOLD <- 0.05

table(genes.expression.small$adj.P.Val.Aph <= FDR_THRESHOLD) 
table(genes.expression.small$adj.P.Val.Aph_KH7 <= FDR_THRESHOLD) 
table(genes.expression.small$adj.P.Val.KH7 <= FDR_THRESHOLD) 
table(genes.expression.small$adj.P.Val.Noc <= FDR_THRESHOLD) 


# genes.expression.small$Selection <- ifelse(( (genes.expression.small$adj.P.Val.Aph < FDR_THRESHOLD) || 
#                                             (genes.expression.small$adj.P.Val.Aph_KH7 < FDR_THRESHOLD) ||
#                                             (genes.expression.small$adj.P.Val.KH7 < FDR_THRESHOLD) ||
#                                             (genes.expression.small$adj.P.Val.Noc < FDR_THRESHOLD) ),
#                                             "TRUE", "FALSE") 

# table(genes.expression.small$Selection)  ##  TRUE 12956 
# length(genes.expression.small$Selection) ## 12956

########################################################################################################################
########################################################################################################################
####### printing the FILE with the SELECTED values :

# write.table(genes.expression.small, 
#            file=paste(NAME, ".file.with.SELECTION.txt", sep=""), 
#            sep="\t") 

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## starting to use 
######################################################################################################################## THE LARGE MATRIX

####### to use two datasets in order to retrieve the genes :
####### genes.expression.large.tpm
####### genes.expression.large.fpkm
####### genes.expression.large

head(genes.expression.large) ### using GENE_NAME
dim(genes.expression.large)  ### using GENE_NAME
tail(genes.expression.large) ### using GENE_NAME

########################################################################################################################
########################################################################################################################
##################################################################################################
##################################################################################################

# > length(genes.expression.large$GENE_NAME)
# [1] 58381
# > length(unique(genes.expression.large$GENE_NAME))
# [1] 56832
#
# X <- genes.expression.large
# Y <- subset(X, !duplicated(X$GENE_NAME))
# dim(Y) ### [1] 56832

genes.expression.large.unique <- subset(genes.expression.large, 
                                        !duplicated(genes.expression.large$GENE_NAME)) 

dim(genes.expression.large.unique)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
##############################################################################  the PCA/MDS analysis on CELL_CYCLE_GENES

genes_cell_cyle <- read.delim("genes.KEGG_Cell_Cycle_GENES.txt", 
                              header=TRUE, sep="\t", stringsAsFactors=F)

genes_cell_cyle$GENE_NAME <- genes_cell_cyle$Gene

head(genes_cell_cyle)
tail(genes_cell_cyle)
dim(genes_cell_cyle)

# genes_cell_cycle_and_info <- merge(genes_cell_cyle, 
#                                  genes.expression.large.unique, 
#                                  by.x = Gene, 
#                                  by.y = GENE_NAME, 
#                                  all.x = TRUE)

genes_cell_cycle_and_info <- join(genes_cell_cyle, genes.expression.large.unique, type = "inner")

head(genes_cell_cycle_and_info)
tail(genes_cell_cycle_and_info)
dim(genes_cell_cycle_and_info)

write.table( genes_cell_cycle_and_info,
             file="genes.KEGG_Cell_Cycle_GENES.with.info.expression.txt", 
             sep="\t", 
             quote = FALSE, 
             row.names = FALSE,
             col.names = TRUE) 
 
##################################################################################################
##################################################################################################

genes_cell_cycle_and_info.fpkm <- subset(genes_cell_cycle_and_info, 
                                     select=c("GENE_NAME", 
                                     "DMSO1_lane1.FPKM", "DMSO1_lane2.FPKM",
                                     "DMSO2_lane1.FPKM", "DMSO2_lane2.FPKM",
                                     "DMSO3_lane1.FPKM", "DMSO3_lane2.FPKM", 
                                     "Aph1.FPKM", "Aph2.FPKM", "Aph3.FPKM", 
                                     "Aph_KH7_1.FPKM","Aph_KH7_2.FPKM","Aph_KH7_3.FPKM",
                                     "KH7_1.FPKM", "KH7_2.FPKM", "KH7_3.FPKM",
                                     "Noc_1.FPKM", "Noc_2.FPKM", "Noc_3.FPKM" ), 
                                      na.rm = TRUE)

rownames(genes_cell_cycle_and_info.fpkm) <- genes_cell_cycle_and_info.fpkm$GENE_NAME 
genes_cell_cycle_and_info.fpkm <- genes_cell_cycle_and_info.fpkm[,-1]

head(genes_cell_cycle_and_info.fpkm)
dim(genes_cell_cycle_and_info.fpkm)

#####################################################################################################
#####################################################################################################
#### R introduces NA

# genes_cell_cycle_and_info.NA <- subset(genes_cell_cycle_and_info, is.na(CHR))
# dim(genes_cell_cycle_and_info.NA)      

# genes_cell_cycle_and_info.non.NA <- subset(genes_cell_cycle_and_info, !is.na(CHR))
# dim(genes_cell_cycle_and_info.non.NA)  

# genes_cell_cycle_and_info <- genes_cell_cycle_and_info.non.NA
# dim(genes_cell_cycle_and_info)

##################################################################################################
##################################################################################################
################################################################################################## HEATMAP

pdf(paste("genes.KEGG_Cell_Cycle_GENES.with.info.expression.heatmap.display.FPKM.pdf", sep=""))
par(las=2)
par(mar=c(8,4,2,2))
par(cex.main=0.6)
    heatmap.2(as.matrix(genes_cell_cycle_and_info.fpkm), col=bluered(149),
                        scale="row",trace="none",
                        cexRow=0.6, cexCol=0.6, cex.main=0.6,
                        Rowv=TRUE, symkey=FALSE, labRow=NA,
                        key=T, keysize=1.5, density.info="none",
                        main="heatmap of KEGG_Cell_Cycle genes")
dev.off()

##################################################################################################
##################################################################################################
##################################################################################################

png(paste("genes.KEGG_Cell_Cycle_GENES.with.info.expression.heatmap.display.FPKM.png", sep=""))
par(las=2)
par(mar=c(8,4,2,2))
par(cex.main=0.6)
    heatmap.2(as.matrix(genes_cell_cycle_and_info.fpkm), col=bluered(149),
                        scale="row",trace="none",
                        cexRow=0.6, cexCol=0.6, cex.main=0.6,
                        Rowv=TRUE, symkey=FALSE, labRow=NA,
                        key=T, keysize=1.5, density.info="none",
                        main="heatmap of KEGG_Cell_Cycle genes")
dev.off()

##################################################################################################
##################################################################################################
################################################################################################## PCA
################################################################################################## the PCA analysis :

group <- factor( c(rep("DMSO",6), rep("Aph", 3), rep("Aph_KH7", 3), 
                   rep("KH7", 3), rep("Noc", 3)) )

pca <- prcomp(t(genes_cell_cycle_and_info.tpm))

###### plotting the pca$x : in PDF format :

color = c(rep("red",6), rep("orange",3), rep("green",3), rep("blue",3), rep("violet",3))
 
pdf(paste("genes.KEGG_Cell_Cycle_GENES.with.info.expression.PCA.display.in.3D.pdf", sep=""), width=10, height=10)

scatterplot3d(pca$x[,1:3],
                           color=colors,
                           pch=18,
                           main="PCA analysis of Cell Cycle Genes", 
                           grid=TRUE, 
                           box=TRUE)

legend(# s3d$xyz.convert(7, 4, 4), 
       "left",
       legend = row.names(pca$x[,1:3]),
       col = colors,                                     ### to use COL, and not COLORS
       pch = 16)

dev.off()

##################################################################################################
##################################################################################################
##################################################################################################
############################################################## plotting the pca$x : in PNG format :

color = c(rep("red",6), rep("orange",3), rep("green",3), rep("blue",3), rep("violet",3))
 
png(paste("genes.KEGG_Cell_Cycle_GENES.with.info.expression.PCA.display.in.3D.png", sep=""))

scatterplot3d(pca$x[,1:3],
                           color=colors,
                           pch=18,
                           main="PCA analysis of Cell Cycle Genes", 
                           grid=TRUE, 
                           box=TRUE)

legend(# s3d$xyz.convert(7, 4, 4), 
       "left",
       legend = row.names(pca$x[,1:3]),
       col = colors,                                     ### to use COL, and not COLORS
       pch = 16)

dev.off()

######################################################################################
######################################################################################
######################################################################################
######################################################################################

pca.df <- data.frame(PCA1=pca$x[,1], 
                     PCA2=pca$x[,2], 
                     PCA3=pca$x[,3], 
                     group=group)

## Here we are plotting PCA1 vs PCA2

pdf(paste("genes.KEGG_Cell_Cycle_GENES.with.info.expression.PCA.display.PC1.vs.PC2.pdf", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

##############################################

png(paste("genes.KEGG_Cell_Cycle_GENES.with.info.expression.PCA.display.PC1.vs.PC2.png", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## to make the BOXPLOTS of log2FC :
# the columns that we are interested in :

# logFC.Aph
# logFC.Aph_KH7
# logFC.KH7
# logFC.Noc

min(genes_cell_cycle_and_info$logFC.Aph, na.rm = TRUE)
max(genes_cell_cycle_and_info$logFC.Aph, na.rm = TRUE)

min(genes_cell_cycle_and_info$logFC.Aph_KH7, na.rm = TRUE)
max(genes_cell_cycle_and_info$logFC.Aph_KH7, na.rm = TRUE)

min(genes_cell_cycle_and_info$logFC.KH7, na.rm = TRUE)
max(genes_cell_cycle_and_info$logFC.KH7, na.rm = TRUE)

min(genes_cell_cycle_and_info$logFC.Noc, na.rm = TRUE)
max(genes_cell_cycle_and_info$logFC.Noc, na.rm = TRUE)

# we can choose a scale from -3.5 to 3.5 on LOG2 scale

# colnames(genes_cell_cycle_and_info)
# [1] "Gene"              "GENE_NAME"         "GENE_NAME_ID"     
# [4] "CHR"               "START"             "END"              
# [7] "STRAND"            "GENE_ID"           "GENE_TYPE"        
#[10] "DMSO1_lane1.count" "DMSO1_lane1.TPM"   "DMSO1_lane1.FPKM" 
#[13] "DMSO1_lane2.count" "DMSO1_lane2.TPM"   "DMSO1_lane2.FPKM" 
#[16] "DMSO2_lane1.count" "DMSO2_lane1.TPM"   "DMSO2_lane1.FPKM" 
#[19] "DMSO2_lane2.count" "DMSO2_lane2.TPM"   "DMSO2_lane2.FPKM" 
#[22] "DMSO3_lane1.count" "DMSO3_lane1.TPM"   "DMSO3_lane1.FPKM" 
#[25] "DMSO3_lane2.count" "DMSO3_lane2.TPM"   "DMSO3_lane2.FPKM" 
#[28] "Aph1.count"        "Aph1.TPM"          "Aph1.FPKM"        
#[31] "Aph2.count"        "Aph2.TPM"          "Aph2.FPKM"        
#[34] "Aph3.count"        "Aph3.TPM"          "Aph3.FPKM"        
#[37] "Aph_KH7_1.count"   "Aph_KH7_1.TPM"     "Aph_KH7_1.FPKM"   
#[40] "Aph_KH7_2.count"   "Aph_KH7_2.TPM"     "Aph_KH7_2.FPKM"   
#[43] "Aph_KH7_3.count"   "Aph_KH7_3.TPM"     "Aph_KH7_3.FPKM"   
#[46] "KH7_1.count"       "KH7_1.TPM"         "KH7_1.FPKM"       
#[49] "KH7_2.count"       "KH7_2.TPM"         "KH7_2.FPKM"       
#[52] "KH7_3.count"       "KH7_3.TPM"         "KH7_3.FPKM"       
#[55] "Noc_1.count"       "Noc_1.TPM"         "Noc_1.FPKM"       
#[58] "Noc_2.count"       "Noc_2.TPM"         "Noc_2.FPKM"       
#[61] "Noc_3.count"       "Noc_3.TPM"         "Noc_3.FPKM"       
#[64] "ID"                "logFC.Aph"         "AveExpr.Aph"      
#[67] "t.Aph"             "P.Value.Aph"       "adj.P.Val.Aph"    
#[70] "B.Aph"             "Gene.Aph"          "ID.Aph"           
#[73] "logFC.Aph_KH7"     "AveExpr.Aph_KH7"   "t.Aph_KH7"        
#[76] "P.Value.Aph_KH7"   "adj.P.Val.Aph_KH7" "B.Aph_KH7"        
#[79] "Gene.Aph_KH7"      "ID.Aph_KH7"        "logFC.KH7"        
#[82] "AveExpr.KH7"       "t.KH7"             "P.Value.KH7"      
#[85] "adj.P.Val.KH7"     "B.KH7"             "Gene.KH7"         
#[88] "ID.KH7"            "logFC.Noc"         "AveExpr.Noc"      
#[91] "t.Noc"             "P.Value.Noc"       "adj.P.Val.Noc"    
#[94] "B.Noc"             "Gene.Noc"          "ID.Noc" 

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## making the BARPLOTS in ggplot2 :

## genes_cell_cycle_and_info$logFC.Aph
## genes_cell_cycle_and_info$logFC.Aph_KH7
## genes_cell_cycle_and_info$logFC.KH7
## genes_cell_cycle_and_info$logFC.Noc


for (i in 1:dim(genes_cell_cycle_and_info)[1])
 {
    # png(paste("cell.cycle.BOXPLOT", genes_cell_cycle_and_info$Gene[i], "the.boxplots.of.log2FC",  "png", sep="."))
    # boxplot( genes_cell_cycle_and_info$logFC.Noc[i],
    #         genes_cell_cycle_and_info$logFC.Aph_KH7[i],
    #         genes_cell_cycle_and_info$logFC.Aph[i],
    #         genes_cell_cycle_and_info$logFC.KH7[i],
    #         ylim=c(-3.5, 3.5), col=c("blue", "red", "green", "brown") )
    # dev.off()

    png(paste("cell.cycle.BARPLOT", genes_cell_cycle_and_info$Gene[i], "the.barplots.of.log2FC",  "png", sep="."), width = 300, height = 600)
    barplot( c(genes_cell_cycle_and_info$logFC.Noc[i],
               genes_cell_cycle_and_info$logFC.Aph_KH7[i],
               genes_cell_cycle_and_info$logFC.Aph[i],
               genes_cell_cycle_and_info$logFC.KH7[i]),
               ylim=c(-3.5, 3.5), col=c("blue", "red", "green", "brown"), 
               # xlab="samples",
               # xlab=c("Noc", "Aph_KH7", "Aph", "KH7"), 
               ylab="log2FC (drug vs control)", 
               main=genes_cell_cycle_and_info$Gene[i])
    axis(side=1, at=1:4, labels=c("Noc", "Aph_KH7", "Aph", "KH7"), las=2)
    dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
############################################################################## in ggplot2, that includes all the SAMPLES 

    df <- data.frame(exp=c("Noc", "Aph_KH7", "Aph", "KH7"), 
                     log2FC=c(genes_cell_cycle_and_info$logFC.Noc[i],
                              genes_cell_cycle_and_info$logFC.Aph_KH7[i],
                              genes_cell_cycle_and_info$logFC.Aph[i],
                              genes_cell_cycle_and_info$logFC.KH7[i]) ) 

    df$exp <- factor(df$exp, levels = df$exp) 

    p <- ggplot(data=df, aes(x=exp, y=log2FC)) +
                geom_bar(stat="identity", fill=c("blue", "red", "green", "brown")) +
                # geom_text(aes(label=exp), vjust=1.6, color="white", size=3.5) +
                # coord_flip() +
                # theme_minimal() +
                # theme(legend.position="top")
                # theme_bw() + 
                # theme_classic() +
                 ggtitle(genes_cell_cycle_and_info$Gene[i]) +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                       text = element_text(size=6), 
                       panel.background = element_rect(fill = "white"))  
   
 
    ggsave(paste("cell.cycle.BARPLOT", genes_cell_cycle_and_info$Gene[i], "the.barplots.of.log2FC.v2",  "png", sep="."), 
           units="cm", width = 3, height = 7)

########################################################################################################################
########################################################################################################################    
########################################################################################################################
############################################################################# in ggplot2, it does not include KH7

    df2 <- data.frame(exp=c("Noc", "Aph_KH7", "Aph"), 
                     log2FC=c(genes_cell_cycle_and_info$logFC.Noc[i],
                              genes_cell_cycle_and_info$logFC.Aph_KH7[i],
                              genes_cell_cycle_and_info$logFC.Aph[i]) )

    df2$exp <- factor(df2$exp, levels = df2$exp) 

    p2 <- ggplot(data=df2, aes(x=exp, y=log2FC)) +
                geom_bar(stat="identity", fill=c("blue", "red", "green")) +
                # geom_text(aes(label=exp), vjust=1.6, color="white", size=3.5) +
                # coord_flip() +
                # theme_minimal() +
                # theme(legend.position="top")
                # theme_bw() + 
                # theme_classic() +
                 ggtitle(genes_cell_cycle_and_info$Gene[i]) +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                       text = element_text(size=6), 
                       panel.background = element_rect(fill = "white"))  
   
 
    ggsave(paste("cell.cycle.BARPLOT", genes_cell_cycle_and_info$Gene[i], "the.barplots.of.log2FC.v3",  "png", sep="."), 
           units="cm", width = 3, height = 7)


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## 
########################################################################################################################
######################################################################################################################## LPS REACTIVE GENES
########################################################################################################################
#################################################################### the PCA/MDS analysis on genes that are LPS_reactive

genes_LPS_reactive <- read.delim("genes.Top50changes_in_LPS_reactive_astrocytes_symbol_HUGO", 
                                  header=T, sep="\t", stringsAsFactors=F)

genes_LPS_reactive$GENE_NAME <- genes_LPS_reactive$Approved_symbol

dim(genes_LPS_reactive)
head(genes_LPS_reactive)
tail(genes_LPS_reactive)

# genes_LPS_reactive_and_info <- merge(genes_LPS_reactive, 
#                                     genes.expression.large.unique, 
#                                     by.x = GENE_NAME, 
#                                     by.y = GENE_NAME, 
#                                     all.x = TRUE)

genes_LPS_reactive_and_info <- join(genes_LPS_reactive, genes.expression.large.unique, type = "inner")

head(genes_LPS_reactive_and_info)
tail(genes_LPS_reactive_and_info)
dim(genes_LPS_reactive_and_info)

write.table( genes_LPS_reactive_and_info,
             file="genes.LPS_reactive.with.info.expression.txt", 
             sep="\t", 
             quote = FALSE, 
             row.names = FALSE,
             col.names = TRUE) 

#####################################################################################################
#####################################################################################################

# genes_LPS_reactive_and_info.NA <- subset(genes_LPS_reactive_and_info, is.na(CHR))
# dim(genes_LPS_reactive_and_info.NA)      ### 

# genes_LPS_reactive_and_info.non.NA <- subset(genes_LPS_reactive_and_info, !is.na(CHR))
# dim(genes_LPS_reactive_and_info.non.NA)  ###

# genes_LPS_reactive_and_info <- genes_LPS_reactive_and_info.non.NA
# dim(genes_LPS_reactive_and_info)

#####################################################################################################
#####################################################################################################

genes_LPS_reactive_and_info.fpkm <-   subset(genes_LPS_reactive_and_info, 
                                      select=c("GENE_NAME", 
                                     "DMSO1_lane1.FPKM", "DMSO1_lane2.FPKM",
                                     "DMSO2_lane1.FPKM", "DMSO2_lane2.FPKM",
                                     "DMSO3_lane1.FPKM", "DMSO3_lane2.FPKM", 
                                     "Aph1.FPKM", "Aph2.FPKM", "Aph3.FPKM", 
                                     "Aph_KH7_1.FPKM","Aph_KH7_2.FPKM","Aph_KH7_3.FPKM",
                                     "KH7_1.FPKM", "KH7_2.FPKM", "KH7_3.FPKM",
                                     "Noc_1.FPKM", "Noc_2.FPKM", "Noc_3.FPKM" ), 
                                      na.rm = TRUE)

rownames(genes_LPS_reactive_and_info.fpkm) <- genes_LPS_reactive_and_info.fpkm$GENE_NAME 
genes_LPS_reactive_and_info.fpkm <- genes_LPS_reactive_and_info.fpkm[,-1]

head(genes_LPS_reactive_and_info.fpkm)
dim(genes_LPS_reactive_and_info.fpkm)

################################################################################################## HEATMAP
################################################################################################## analysis
##################################################################################################

pdf(paste("genes.LPS_reactive.with.info.expression.heatmap.display.FPKM.pdf", sep=""))
par(las=2)
par(mar=c(8,4,2,2))
par(cex.main=0.6)
    heatmap.2(as.matrix(genes_LPS_reactive_and_info.fpkm), col=bluered(149),
                        scale="row",trace="none",
                        cexRow=0.6, cexCol=0.6, cex.main=0.6,
                        Rowv=TRUE, symkey=FALSE, labRow=NA,
                        key=T, keysize=1.5, density.info="none",
                        main="heatmap of LPS reactive genes")
dev.off()

##################################################################################################
##################################################################################################

png(paste("genes.LPS_reactive.with.info.expression.heatmap.display.FPKM.png", sep=""))
par(las=2)
par(mar=c(8,4,2,2))
par(cex.main=0.6)
    heatmap.2(as.matrix(genes_LPS_reactive_and_info.fpkm), col=bluered(149),
                        scale="row",trace="none",
                        cexRow=0.6, cexCol=0.6, cex.main=0.6,
                        Rowv=TRUE, symkey=FALSE, labRow=NA,
                        key=T, keysize=1.5, density.info="none",
                        main="heatmap of LPS reactive genes")
dev.off()

##################################################################################################
################################################################################################## PCA
################################################################################################## analysis

group <- factor( c(rep("DMSO",6), rep("Aph", 3), rep("Aph_KH7", 3), 
                   rep("KH7", 3), rep("Noc", 3)) )

pca <- prcomp(t(genes_LPS_reactive_and_info.tpm))

colors = c(rep("red",6), rep("orange",3), rep("green",3), rep("blue",3), rep("violet",3))

#########################################################################################
######################################################################################### plotting the pca$x :
 
pdf(paste("genes.LPS_reactive.with.info.expression.PCA.display.in.3D.pdf", sep=""))
    scatterplot3d(pca$x[,1:3],
                           color=colors,
                           pch=18,
                           main="PCA analysis of LPS reactive genes", 
                           grid=TRUE, 
                           box=TRUE)

    legend(# s3d$xyz.convert(7, 4, 4), 
           "left",
           legend = row.names(pca$x[,1:3]),
           col = colors,                                     ### it has to use COL, and not COLORS
           pch = 16)

dev.off()

#########################################################################################
######################################################################################### plotting the pca$x :
 
png(paste("genes.LPS_reactive.with.info.expression.PCA.display.in.3D.png", sep=""))
    scatterplot3d(pca$x[,1:3],
                           color=colors,
                           pch=18,
                           main="PCA analysis of LPS reactive genes", 
                           grid=TRUE, 
                           box=TRUE)

    legend(# s3d$xyz.convert(7, 4, 4), 
           "left",
           legend = row.names(pca$x[,1:3]),
           col = colors,                                     ### it has to use COL, and not COLORS
           pch = 16)

dev.off()

######################################################################################
###################################################################################### PCA1 and PCA2

pca.df <- data.frame(PCA1=pca$x[,1], 
                     PCA2=pca$x[,2], 
                     PCA3=pca$x[,3], 
                     group=group)

######################################################################################
############################################################################## plotting PCA1 vs PCA2
######################################################################################

pdf(paste("genes.LPS_reactive.with.info.expression.PCA.display.PC1.vs.PC2.pdf", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

png(paste("genes.LPS_reactive.with.info.expression.PCA.display.PC1.vs.PC2.png", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## the BARPLOTS in ggplot2 :

## genes_LPS_reactive_and_info$logFC.Aph
## genes_LPS_reactive_and_info$logFC.Aph_KH7
## genes_LPS_reactive_and_info$logFC.KH7
## genes_LPS_reactive_and_info$logFC.Noc

dim(genes_LPS_reactive_and_info)[1]

for (i in 1:dim(genes_LPS_reactive_and_info)[1])
 {
    png(paste("LPS.reactive.BARPLOT", genes_LPS_reactive_and_info$GENE_NAME[i], "the.barplots.of.log2FC",  "png", sep="."), width = 300, height = 600)
    barplot( c(genes_LPS_reactive_and_info$logFC.Noc[i],
               genes_LPS_reactive_and_info$logFC.Aph_KH7[i],
               genes_LPS_reactive_and_info$logFC.Aph[i],
               genes_LPS_reactive_and_info$logFC.KH7[i]),
               ylim=c(-3.5, 3.5), col=c("blue", "red", "green", "brown"), 
               # xlab="samples",
               # xlab=c("Noc", "Aph_KH7", "Aph", "KH7"), 
               ylab="log2FC (drug vs control)", 
               main=genes_LPS_reactive_and_info$GENE_NAME[i])
    axis(side=1, at=1:4, labels=c("Noc", "Aph_KH7", "Aph", "KH7"), las=2)
    dev.off()

######################################################################################################################
################################################################ the version in ggplot2, that includes all the SAMPLES 

    df <- data.frame(exp=c("Noc", "Aph_KH7", "Aph", "KH7"), 
                     log2FC=c(genes_LPS_reactive_and_info$logFC.Noc[i],
                              genes_LPS_reactive_and_info$logFC.Aph_KH7[i],
                              genes_LPS_reactive_and_info$logFC.Aph[i],
                              genes_LPS_reactive_and_info$logFC.KH7[i]) ) 

    df$exp <- factor(df$exp, levels = df$exp) 

    p <- ggplot(data=df, aes(x=exp, y=log2FC)) +
                geom_bar(stat="identity", fill=c("blue", "red", "green", "brown")) +
                # geom_text(aes(label=exp), vjust=1.6, color="white", size=3.5) +
                # coord_flip() +
                # theme_minimal() +
                # theme(legend.position="top")
                # theme_bw() + 
                # theme_classic() +
                 ggtitle(genes_LPS_reactive_and_info$GENE_NAME[i]) +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                       text = element_text(size=6), 
                       panel.background = element_rect(fill = "white"))  
   
 
    ggsave(paste("LPS.reactive.BARPLOT", genes_LPS_reactive_and_info$GENE_NAME[i], "the.barplots.of.log2FC.v2",  "png", sep="."), 
           units="cm", width = 3, height = 7)
    
######################################################################################################################
#################################################################### the version in ggplot2, that does not include KH7

    df2 <- data.frame(exp=c("Noc", "Aph_KH7", "Aph"), 
                     log2FC=c(genes_LPS_reactive_and_info$logFC.Noc[i],
                              genes_LPS_reactive_and_info$logFC.Aph_KH7[i],
                              genes_LPS_reactive_and_info$logFC.Aph[i]) )

    df2$exp <- factor(df2$exp, levels = df2$exp) 

    p2 <- ggplot(data=df2, aes(x=exp, y=log2FC)) +
                geom_bar(stat="identity", fill=c("blue", "red", "green")) +
                # geom_text(aes(label=exp), vjust=1.6, color="white", size=3.5) +
                # coord_flip() +
                # theme_minimal() +
                # theme(legend.position="top")
                # theme_bw() + 
                # theme_classic() +
                 ggtitle(genes_LPS_reactive_and_info$GENE_NAME[i]) +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                       text = element_text(size=6), 
                       panel.background = element_rect(fill = "white"))  
   
 
    ggsave(paste("LPS.reactive.BARPLOT", genes_LPS_reactive_and_info$GENE_NAME[i], "the.barplots.of.log2FC.v3",  "png", sep="."), 
           units="cm", width = 3, height = 7)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## MCAO REACTIVE genes
########################################################################################################################
################################################################################## doing the PCA/MDS analysis on genes that are MCAO_reactive

genes_MCAO_reactive <- read.delim("genes.Top50changes_in_MCAO_reactive_astrocytes_symbol_HUGO_v2", 
                                  header=T, sep="\t", stringsAsFactors=F)

genes_MCAO_reactive$GENE_NAME <- genes_MCAO_reactive$Approved_symbol

dim(genes_MCAO_reactive)
head(genes_MCAO_reactive)
tail(genes_MCAO_reactive)

########################################################################################################################
########################################################################################################################
########################################################################################################################

genes_MCAO_reactive_and_info <- join(genes_MCAO_reactive, genes.expression.large.unique, type = "inner")

dim(genes_MCAO_reactive_and_info)
head(genes_MCAO_reactive_and_info)
tail(genes_MCAO_reactive_and_info)

########################################################################################################################
########################################################################################################################
########################################################################################################################


write.table( genes_MCAO_reactive_and_info,
             file="genes.MCAO_reactive.with.info.expression.txt", 
             sep="\t", 
             quote = FALSE, 
             row.names = FALSE,
             col.names = TRUE) 

#####################################################################################################

# genes_MCAO_reactive_and_info.NA <- subset(genes_MCAO_reactive_and_info, is.na(CHR))
# dim(genes_MCAO_reactive_and_info.NA)      ### 

# genes_MCAO_reactive_and_info.non.NA <- subset(genes_MCAO_reactive_and_info, !is.na(CHR))
# dim(genes_MCAO_reactive_and_info.non.NA)  ###

# genes_MCAO_reactive_and_info <- genes_MCAO_reactive_and_info.non.NA
# dim(genes_MCAO_reactive_and_info)

######################################################################################################

genes_MCAO_reactive_and_info.fpkm <-   subset(genes_MCAO_reactive_and_info, 
                                             select=c("GENE_NAME", 
                                            "DMSO1_lane1.FPKM", "DMSO1_lane2.FPKM",
                                            "DMSO2_lane1.FPKM", "DMSO2_lane2.FPKM",
                                            "DMSO3_lane1.FPKM", "DMSO3_lane2.FPKM", 
                                            "Aph1.FPKM", "Aph2.FPKM", "Aph3.FPKM", 
                                            "Aph_KH7_1.FPKM","Aph_KH7_2.FPKM","Aph_KH7_3.FPKM",
                                            "KH7_1.FPKM", "KH7_2.FPKM", "KH7_3.FPKM",
                                            "Noc_1.FPKM", "Noc_2.FPKM", "Noc_3.FPKM" ), 
                                             na.rm = TRUE)

rownames(genes_MCAO_reactive_and_info.fpkm) <- genes_MCAO_reactive_and_info.fpkm$GENE_NAME
genes_MCAO_reactive_and_info.fpkm <- genes_MCAO_reactive_and_info.fpkm[,-1]

head(genes_MCAO_reactive_and_info.fpkm)
dim(genes_MCAO_reactive_and_info.fpkm)

##################################################################################################
##################################################################################################
##################################################################################################
################################################################################################## the HEATMAP analysis :

pdf(paste("genes.MCAO_reactive.with.info.expression.heatmap.display.FPKM.pdf", sep=""))
par(las=2)
par(mar=c(8,4,2,2))
par(cex.main=0.6)
    heatmap.2(as.matrix(genes_MCAO_reactive_and_info.fpkm), col=bluered(149),
                        scale="row",trace="none",
                        cexRow=0.6, cexCol=0.6, cex.main=0.6,
                        Rowv=TRUE, symkey=FALSE, labRow=NA,
                        key=T, keysize=1.5, density.info="none",
                        main="heatmap of MCAO reactive genes")
dev.off()

###################################################################################################

png(paste("genes.MCAO_reactive.with.info.expression.heatmap.display.FPKM.png", sep=""))
par(las=2)
par(mar=c(8,4,2,2))
par(cex.main=0.6)
    heatmap.2(as.matrix(genes_MCAO_reactive_and_info.fpkm), col=bluered(149),
                        scale="row",trace="none",
                        cexRow=0.6, cexCol=0.6, cex.main=0.6,
                        Rowv=TRUE, symkey=FALSE, labRow=NA,
                        key=T, keysize=1.5, density.info="none",
                        main="heatmap of MCAO reactive genes")
dev.off()

##################################################################################################
##################################################################################################
################################################################################################## PCA

group <- factor( c(rep("DMSO",6), rep("Aph", 3), rep("Aph_KH7", 3), rep("KH7", 3), rep("Noc", 3)) )

pca <- prcomp(t(genes_LPS_reactive_and_info.tpm))

colors = c(rep("red",6), rep("orange",3), rep("green",3), rep("blue",3), rep("violet",3))

################################################################################# plotting the pca$x :
 
pdf(paste("genes.MCAO_reactive.with.info.expression.PCA.display.in.3D.pdf", sep=""))
s3d <- scatterplot3d(pca$x[,1:3],
                           color = colors, 
                           pch=18,
                           main="PCA analysis of MCAO reactive genes", 
                           grid=TRUE, 
                           box=TRUE)

legend(# s3d$xyz.convert(7, 4, 4), 
       "left",
       legend = row.names(pca$x[,1:3]),
       col = colors,                                     ### it has to use COL, and not COLORS
       pch = 16)

dev.off()

################################################################################## plotting the pca$x :
 
png(paste("genes.MCAO_reactive.with.info.expression.PCA.display.in.3D.png", sep=""))
s3d <- scatterplot3d(pca$x[,1:3],
                           color = colors, 
                           pch=18,
                           main="PCA analysis of MCAO reactive genes", 
                           grid=TRUE, 
                           box=TRUE)

legend(# s3d$xyz.convert(7, 4, 4), 
       "left",
       legend = row.names(pca$x[,1:3]),
       col = colors,                                     ### it has to use COL, and not COLORS
       pch = 16)

dev.off()

######################################################################################
######################################################################################

pca.df <- data.frame(PCA1=pca$x[,1], 
                     PCA2=pca$x[,2], 
                     PCA3=pca$x[,3], 
                     group=group)

###################################################################################### plotting PCA1 vs PCA2

pdf(paste("genes.MCAO_reactive.with.info.expression.PCA.display.PC1.vs.PC2.pdf", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

png(paste("genes.MCAO_reactive.with.info.expression.PCA.display.PC1.vs.PC2.png", sep=""))
ggplot(pca.df, 
       aes(x=PCA1, y=PCA2, color=group, label=rownames(pca.df))) +
       geom_point(size=3) +
       # geom_text(col='black', size=4) +
       theme_bw() +
       theme(legend.position="top", 
             legend.title=element_blank(), 
             legend.key = element_blank()) +
             labs(x="PC1", y="PC2") + 
        ggtitle("PCA analysis : PC1 vs PC2")
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## the BARPLOTS in ggplot2 :

## genes_MCAO_reactive_and_info$logFC.Aph
## genes_MCAO_reactive_and_info$logFC.Aph_KH7
## genes_MCAO_reactive_and_info$logFC.KH7
## genes_MCAO_reactive_and_info$logFC.Noc

dim(genes_MCAO_reactive_and_info)[1]

for (i in 1:dim(genes_MCAO_reactive_and_info)[1])
 {

    png(paste("MCAO.reactive.BARPLOT", genes_MCAO_reactive_and_info$GENE_NAME[i], "the.barplots.of.log2FC",  "png", sep="."), width = 300, height = 600)
    barplot( c(genes_MCAO_reactive_and_info$logFC.Noc[i],
               genes_MCAO_reactive_and_info$logFC.Aph_KH7[i],
               genes_MCAO_reactive_and_info$logFC.Aph[i],
               genes_MCAO_reactive_and_info$logFC.KH7[i]),
               ylim=c(-3.5, 3.5), col=c("blue", "red", "green", "brown"), 
               # xlab="samples",
               # xlab=c("Noc", "Aph_KH7", "Aph", "KH7"), 
               ylab="log2FC (drug vs control)", 
               main=genes_MCAO_reactive_and_info$GENE_NAME[i])
    axis(side=1, at=1:4, labels=c("Noc", "Aph_KH7", "Aph", "KH7"), las=2)
    dev.off()

#######################################################################################################################
############################################################################################# the version in ggplot2, that includes all the SAMPLES 

    df <- data.frame(exp=c("Noc", "Aph_KH7", "Aph", "KH7"), 
                     log2FC=c(genes_MCAO_reactive_and_info$logFC.Noc[i],
                              genes_MCAO_reactive_and_info$logFC.Aph_KH7[i],
                              genes_MCAO_reactive_and_info$logFC.Aph[i],
                              genes_MCAO_reactive_and_info$logFC.KH7[i]) ) 

    df$exp <- factor(df$exp, levels = df$exp) 

    p <- ggplot(data=df, aes(x=exp, y=log2FC)) +
                geom_bar(stat="identity", fill=c("blue", "red", "green", "brown")) +
                # geom_text(aes(label=exp), vjust=1.6, color="white", size=3.5) +
                # coord_flip() +
                # theme_minimal() +
                # theme(legend.position="top")
                # theme_bw() + 
                # theme_classic() +
                 ggtitle(genes_MCAO_reactive_and_info$GENE_NAME[i]) +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                       text = element_text(size=6), 
                       panel.background = element_rect(fill = "white"))  
   
 
    ggsave(paste("MCAO.reactive.BARPLOT", genes_MCAO_reactive_and_info$GENE_NAME[i], "the.barplots.of.log2FC.v2",  "png", sep="."), 
           units="cm", width = 3, height = 7)
           
#######################################################################################################################    
#######################################################################################################################
#################################################################### the version in ggplot2, that doees not include KH7

    df2 <- data.frame(exp=c("Noc", "Aph_KH7", "Aph"), 
                     log2FC=c(genes_MCAO_reactive_and_info$logFC.Noc[i],
                              genes_MCAO_reactive_and_info$logFC.Aph_KH7[i],
                              genes_MCAO_reactive_and_info$logFC.Aph[i]) )

    df2$exp <- factor(df2$exp, levels = df2$exp) 

    p2 <- ggplot(data=df2, aes(x=exp, y=log2FC)) +
                geom_bar(stat="identity", fill=c("blue", "red", "green")) +
                # geom_text(aes(label=exp), vjust=1.6, color="white", size=3.5) +
                # coord_flip() +
                # theme_minimal() +
                # theme(legend.position="top")
                # theme_bw() + 
                # theme_classic() +
                 ggtitle(genes_MCAO_reactive_and_info$GENE_NAME[i]) +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                       text = element_text(size=6), 
                       panel.background = element_rect(fill = "white"))  
   
 
    ggsave(paste("MCAO.reactive.BARPLOT", genes_MCAO_reactive_and_info$GENE_NAME[i], "the.barplots.of.log2FC.v3",  "png", sep="."), 
           units="cm", width = 3, height = 7)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
######################################################################################################################## 
########################################################################################################################
########################################################################################################################
########################################################################################################################
###### the MA PLOTS : a separate file
########################################################################################################################
########################################################################################################################
########################################################################################################################
###### the SCATTER PLOTS : separate file
########################################################################################################################
########################################################################################################################
########################################################################################################################
###### the VOLCANO PLOTS : separate file
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

# an example from : http://www.sthda.com/english/wiki/scatterplot3d-3d-graphics-r-software-and-data-visualization

# http://www.sthda.com/english/wiki/scatterplot3d-3d-graphics-r-software-and-data-visualization

########################################################################################################################
########################################################################################################################
########################################################################################################################
