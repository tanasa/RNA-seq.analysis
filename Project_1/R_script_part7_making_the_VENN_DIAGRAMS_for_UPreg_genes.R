###############################################################################################
###############################################################################################

library("ggplot2")
library("reshape2")
library("data.table")
library("dplyr")
library("plyr")
library("tidyr")
library("gplots")
library("pheatmap")
library("RColorBrewer")
library("scatterplot3d")
library("limma")
library("edgeR")
library("Glimma")
library("DESeq2")

library("VennDiagram") 
library("Vennerable")
library("gplots")

################################################################################################
################################################################################################ UP
# reading the files with UP-regulated genes
# analysis.LIMMA.genes.counts.Aph_KH7.only.DEG.and.UP.info.all.samples
# analysis.LIMMA.genes.counts.Aph.only.DEG.and.UP.info.all.samples
# analysis.LIMMA.genes.counts.KH7.only.DEG.and.UP.info.all.samples
# analysis.LIMMA.genes.counts.Noc.only.DEG.and.UP.info.all.samples

up.Aph <- read.delim("analysis.LIMMA.genes.counts.Aph.only.DEG.and.UP.info.all.samples", header=T, sep="\t", stringsAsFactors=F)
dim(up.Aph)

up.Aph.KH7 <- read.delim("analysis.LIMMA.genes.counts.Aph_KH7.only.DEG.and.UP.info.all.samples", header=T, sep="\t", stringsAsFactors=F)
dim(up.Aph.KH7)

up.KH7 <- read.delim("analysis.LIMMA.genes.counts.KH7.only.DEG.and.UP.info.all.samples", header=T, sep="\t", stringsAsFactors=F)
dim(up.KH7)

up.Noc <- read.delim("analysis.LIMMA.genes.counts.Noc.only.DEG.and.UP.info.all.samples", header=T, sep="\t", stringsAsFactors=F) 
dim(up.Noc)

################################################################################################ about the number of genes :

dim(up.Aph.KH7)
dim(up.KH7)
dim(up.Aph)
dim(up.Noc)

# dim(up.Aph.KH7)
#[1] 3387  103
# dim(up.KH7)
#[1] 3455  103
# dim(up.Aph)
#[1] 2259  103
# dim(up.Noc)
#[1] 3357  103

sets <- list( genes.up.Aph = up.Aph$Gene,
              genes.up.Aph.KH7 = up.Aph.KH7$Gene,
              genes.up.KH7 = up.KH7$Gene, 
              genes.up.Noc = up.Noc$Gene)  

sets_venn <- venn.diagram(sets, filename=NULL)

png(paste("intersect.4sets.VENN.up.reg.png", sep=""), width = 680, height = 680)
grid.newpage()
grid.draw(sets_venn)
dev.off()

################################################################################################
################################################################################################

sets1 <- list( genes.up.Aph.KH7 = up.Aph.KH7$Gene,
               genes.up.KH7 = up.KH7$Gene )  

sets1_venn <- venn.diagram(sets1, filename=NULL)

png(paste("intersect.2sets.VENN.up.reg.Aph_KH7.vs.KH7.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets1_venn)
dev.off()

################################################################################################
################################################################################################

sets2 <- list( genes.up.Aph.KH7 = up.Aph.KH7$Gene,
               genes.up.Aph = up.Aph$Gene )  

sets2_venn <- venn.diagram(sets2, filename=NULL)

png(paste("intersect.2sets.VENN.up.reg.Aph_KH7.vs.Aph.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets2_venn)
dev.off()

################################################################################################
################################################################################################

sets3 <- list( genes.up.Aph.KH7 = up.Aph.KH7$Gene,
               genes.up.Noc = up.Noc$Gene )  

sets3_venn <- venn.diagram(sets3, filename=NULL)

png(paste("intersect.2sets.VENN.up.reg.Aph_KH7.vs.Noc.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets3_venn)
dev.off()

################################################################################################
################################################################################################ intersecting 3 sets

sets4 <- list( genes.up.Aph = up.Aph$Gene,
               genes.up.Aph.KH7 = up.Aph.KH7$Gene,
               genes.up.Noc = up.Noc$Gene)  

sets4_venn <- venn.diagram(sets4, filename=NULL)

png(paste("intersect.3sets.VENN.up.reg.png", sep=""), width = 680, height = 680)
grid.newpage()
grid.draw(sets4_venn)
dev.off()

################################################################################################
################################################################################################

sets5 <- list( genes.up.KH7 = up.KH7$Gene,
               genes.up.Aph = up.Aph$Gene )  

sets5_venn <- venn.diagram(sets5, filename=NULL)

png(paste("intersect.2sets.VENN.up.reg.KH7.vs.Aph.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets5_venn)
dev.off()

################################################################################################
################################################################################################

sets6 <- list( genes.up.KH7 = up.KH7$Gene,
               genes.up.Noc = up.Noc$Gene )  

sets6_venn <- venn.diagram(sets6, filename=NULL)

png(paste("intersect.2sets.VENN.up.reg.KH7.vs.Noc.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets6_venn)
dev.off()

################################################################################################
################################################################################################
################################################################################################
################################################################################################
