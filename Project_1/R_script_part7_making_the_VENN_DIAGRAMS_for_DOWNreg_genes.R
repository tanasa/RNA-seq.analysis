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
################################################################################################ DOWN
# reading the files with DOWN-regulated genes
# analysis.LIMMA.genes.counts.Aph_KH7.only.DEG.and.DOWN.info.all.samples
# analysis.LIMMA.genes.counts.Aph.only.DEG.and.DOWN.info.all.samples
# analysis.LIMMA.genes.counts.KH7.only.DEG.and.DOWN.info.all.samples
# analysis.LIMMA.genes.counts.Noc.only.DEG.and.DOWN.info.all.samples

down.Aph <- read.delim("analysis.LIMMA.genes.counts.Aph.only.DEG.and.DOWN.info.all.samples", header=T, sep="\t", stringsAsFactors=F)
dim(down.Aph)

down.Aph.KH7 <- read.delim("analysis.LIMMA.genes.counts.Aph_KH7.only.DEG.and.DOWN.info.all.samples", header=T, sep="\t", stringsAsFactors=F)
dim(down.Aph.KH7)

down.KH7 <- read.delim("analysis.LIMMA.genes.counts.KH7.only.DEG.and.DOWN.info.all.samples", header=T, sep="\t", stringsAsFactors=F)
dim(down.KH7)

down.Noc <- read.delim("analysis.LIMMA.genes.counts.Noc.only.DEG.and.DOWN.info.all.samples", header=T, sep="\t", stringsAsFactors=F) 
dim(down.Noc)

################################################################################################ about the number of genes :

dim(down.Aph.KH7)
dim(down.KH7)
dim(down.Aph)
dim(down.Noc)

# dim(down.Aph)
#[1] 1844  103
# dim(down.Aph.KH7)
#[1] 3677  103
# dim(down.KH7)
#[1] 3428  103
# dim(down.Noc)
#[1] 3542  10

sets <- list( genes.down.Aph = down.Aph$Gene,
              genes.down.Aph.KH7 = down.Aph.KH7$Gene,
              genes.down.KH7 = down.KH7$Gene, 
              genes.down.Noc = down.Noc$Gene)  

sets_venn <- venn.diagram(sets, filename=NULL)

png(paste("intersect.4sets.VENN.DOWN.reg.png", sep=""), width = 680, height = 680)
grid.newpage()
grid.draw(sets_venn)
dev.off()

################################################################################################
################################################################################################

sets1 <- list( genes.down.Aph.KH7 = down.Aph.KH7$Gene,
               genes.down.KH7 = down.KH7$Gene )  

sets1_venn <- venn.diagram(sets1, filename=NULL)

png(paste("intersect.2sets.VENN.DOWN.reg.Aph_KH7.vs.KH7.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets1_venn)
dev.off()

################################################################################################
################################################################################################

sets2 <- list( genes.down.Aph.KH7 = down.Aph.KH7$Gene,
               genes.down.Aph = down.Aph$Gene )  

sets2_venn <- venn.diagram(sets2, filename=NULL)

png(paste("intersect.2sets.VENN.DOWN.reg.Aph_KH7.vs.Aph.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets2_venn)
dev.off()

################################################################################################
################################################################################################

sets3 <- list( genes.down.Aph.KH7 = down.Aph.KH7$Gene,
               genes.down.Noc = down.Noc$Gene )  

sets3_venn <- venn.diagram(sets3, filename=NULL)

png(paste("intersect.2sets.VENN.DOWN.reg.Aph_KH7.vs.Noc.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets3_venn)
dev.off()

################################################################################################
################################################################################################ intersecting 3 sets

sets4 <- list( genes.down.Aph = down.Aph$Gene,
               genes.down.Aph.KH7 = down.Aph.KH7$Gene,
               genes.down.Noc = down.Noc$Gene)  

sets4_venn <- venn.diagram(sets4, filename=NULL)

png(paste("intersect.3sets.VENN.DOWN.reg.png", sep=""), width = 680, height = 680)
grid.newpage()
grid.draw(sets4_venn)
dev.off()

################################################################################################
################################################################################################

sets5 <- list( genes.down.KH7 = down.KH7$Gene,
               genes.down.Aph = down.Aph$Gene )  

sets5_venn <- venn.diagram(sets5, filename=NULL)

png(paste("intersect.2sets.VENN.DOWN.reg.KH7.vs.Aph.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets5_venn)
dev.off()

################################################################################################
################################################################################################

sets6 <- list( genes.down.KH7 = down.KH7$Gene,
               genes.down.Noc = down.Noc$Gene )  

sets6_venn <- venn.diagram(sets6, filename=NULL)

png(paste("intersect.2sets.VENN.DOWN.reg.KH7.vs.Noc.png", sep=""), width = 680, height = 580)
grid.newpage()
grid.draw(sets6_venn)
dev.off()


################################################################################################
################################################################################################
################################################################################################
################################################################################################
