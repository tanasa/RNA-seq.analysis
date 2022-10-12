######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

library("limma")
library("ggplot2")

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# We are reading these 2 files :
# [4] "RESULTS.limma.compare.CTRL.SOX11MUT.genes.DOWN.forR.txt"   
# [6] "RESULTS.limma.compare.CTRL.SOX11MUT.genes.UP.only.forR.txt"

DOWN <- read.delim("RESULTS.limma.compare.CTRL.SOX11MUT.genes.DOWN.forR.txt" , 
                   header=T, sep="\t", stringsAsFactors=F)

head(DOWN)
dim(DOWN)

UP <- read.delim("RESULTS.limma.compare.CTRL.SOX11MUT.genes.UP.only.forR.txt", 
                   header=T, sep="\t", stringsAsFactors=F)
head(UP)
dim(UP)

###################################################################
###################################################################
### DOWN : we are selecting the genes that are more DOWN-REG by MUT
### [27] "logFC.Control.Sox11"        
### [33] "logFC.Control.Sox11MUT"

DOWN.moreMUT <- subset(DOWN, logFC.Control.Sox11MUT < logFC.Control.Sox11)

# dim(DOWN.moreMUT)

write.table(DOWN.moreMUT, 
            file="genes.DOWN.moreMUT.txt", sep="\t", quote=FALSE, 
            row.names = FALSE, col.names = TRUE)

# in order to show a SCATTER PLOT or BOX PLOT, we generate simpler datasets :

DOWN.moreMUT.FPKM <- subset(DOWN.moreMUT, 
                     select=c("gene_name", "Control_1.FPKM", "Control_2.FPKM", 
                                            "Sox11_1.FPKM", "Sox11_2.FPKM", 
                                            "Sox11MUT_1.FPKM", "Sox11MUT_2.FPKM")) 

png("the.boxplots.of.FPKM.for.DOWN.REG.genes.more.in.MUT.png")
boxplot(DOWN.moreMUT.FPKM$Control_1.FPKM, 
        DOWN.moreMUT.FPKM$Control_2.FPKM, 
        DOWN.moreMUT.FPKM$Sox11_1.FPKM, 
        DOWN.moreMUT.FPKM$Sox11_2.FPKM, 
        DOWN.moreMUT.FPKM$Sox11MUT_1.FPKM,
        DOWN.moreMUT.FPKM$Sox11MUT_2.FPKM, 
        ylim=c(0,50), col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

DOWN.moreMUT.FC <- subset(DOWN.moreMUT, 
                    select=c("gene_name", "logFC.Control.Sox11" , "logFC.Control.Sox11MUT"))

png("the.boxplots.of.FC.for.DOWN.REG.genes.more.in.MUT.png",  width = 200, height = 600)
boxplot(DOWN.moreMUT.FC$logFC.Control.Sox11, 
        DOWN.moreMUT.FC$logFC.Control.Sox11MUT, 
        ylim=c(-5,0), 
        col=c("red", "green"), notch=TRUE, ylab="log2FC", main = "MUT : more DOWN-reg")
dev.off()

wilcox.test(DOWN.moreMUT.FC$logFC.Control.Sox11, DOWN.moreMUT.FC$logFC.Control.Sox11MUT)

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
################################################################## in order to make the SCATTER PLOTS : or to make VOLCANO PLOTS ...
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# or just to select the columns that we are interested in :
# colnames(DOWN.moreMUT) 
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

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###################################################################
###################################################################
### UP : we are selecting the genes that are more UP-REG by MUT
### [27] "logFC.Control.Sox11"        
### [33] "logFC.Control.Sox11MUT"

UP.moreMUT <- subset(UP, logFC.Control.Sox11MUT > logFC.Control.Sox11)

dim(UP.moreMUT)
#[1] 181  38

write.table(UP.moreMUT, 
            file="genes.UP.moreMUT.txt", sep="\t", quote=FALSE, 
            row.names = FALSE, col.names = TRUE)

# in order to show a SCATTER PLOT or a BOX PLOT, we make simpler datasets :

UP.moreMUT.FPKM <- subset(UP.moreMUT, 
                     select=c("gene_name", "Control_1.FPKM", "Control_2.FPKM", 
                                            "Sox11_1.FPKM", "Sox11_2.FPKM", 
                                            "Sox11MUT_1.FPKM", "Sox11MUT_2.FPKM")) 

png("the.boxplots.of.FPKM.for.UP.REG.genes.more.in.MUT.png")
boxplot(UP.moreMUT.FPKM$Control_1.FPKM, 
        UP.moreMUT.FPKM$Control_2.FPKM, 
        UP.moreMUT.FPKM$Sox11_1.FPKM, 
        UP.moreMUT.FPKM$Sox11_2.FPKM, 
        UP.moreMUT.FPKM$Sox11MUT_1.FPKM,
        UP.moreMUT.FPKM$Sox11MUT_2.FPKM,
        ylim=c(0,50), col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

UP.moreMUT.FC <- subset(UP.moreMUT, 
                    select=c("gene_name", "logFC.Control.Sox11" , "logFC.Control.Sox11MUT"))

png("the.boxplots.of.FC.for.UP.REG.genes.more.in.MUT.png",  width = 200, height = 600)
boxplot(UP.moreMUT.FC$logFC.Control.Sox11, 
        UP.moreMUT.FC$logFC.Control.Sox11MUT, 
        ylim=c(0,5), 
        col=c("red", "green"), notch=TRUE, ylab="log2FC", main = "MUT : more UP-reg")
dev.off()

wilcox.test(UP.moreMUT.FC$logFC.Control.Sox11, UP.moreMUT.FC$logFC.Control.Sox11MUT)

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# combining BOTH LISTS :

REG <- rbind(DOWN, UP)

# dim(REG)

# the genes that were selected :

kunche <- read.delim("the_list_of_genes_from_KUNCHE.txt", header=T, sep="\t", stringsAsFactors=F)

# combing the FILES and printing :

kunche.and.reg <- merge(kunche, REG, 
                        by.x = "gene",
                        by.y = "gene_name", 
                        all = FALSE) 

write.table(kunche.and.reg, 
           file="the_list_of_genes_from_KUNCHE.with.EXPRESSION.VALUES.txt",
           sep="\t", quote=FALSE, 
           row.names=FALSE, col.names=TRUE) 

##################################################################################################
##################################################################################################
# in order to display the BOXPLOTS of FC and FPKM : 
##################################################################################################
##################################################################################################
################ we can "activate" the R code below, if/when needed ...

# for (i in 1:dim(kunche.and.reg)[1])
# {
#    png(paste(kunche.and.reg$gene[i], "the.boxplots.of.FPKM",  "png", sep="."))
#    boxplot(kunche.and.reg$Control_1.FPKM[i], 
#            kunche.and.reg$Control_2.FPKM[i], 
#            kunche.and.reg$Sox11_1.FPKM[i], 
#            kunche.and.reg$Sox11_2.FPKM[i], 
#            kunche.and.reg$Sox11MUT_1.FPKM[i],
#            kunche.and.reg$Sox11MUT_2.FPKM[i],
#            ylim=c(0,50), col=c("blue", "blue", "red", "red", "green", "green") )
#   dev.off()

#   if (kunche.and.reg$REGULATION[i] == "UP")
#   {
#      png(paste(kunche.and.reg$gene[i], "the.boxplots.of.log2FC.UP.gene", "png", sep="."), width = 200, height = 600)
#      boxplot(kunche.and.reg$logFC.Control.Sox11[i], 
#              kunche.and.reg$logFC.Control.Sox11MUT[i], 
#              ylim=c(0, 6.5), 
#              col=c("red", "green"), notch=TRUE, ylab="log2FC", main = "")
#      dev.off()
#   }

#   if (kunche.and.reg$REGULATION[i] == "DOWN")
#   {
#      png(paste(kunche.and.reg$gene[i], "the.boxplots.of.log2FC.DOWN.gene", "png", sep="."), width = 200, height = 600)
#      boxplot(kunche.and.reg$logFC.Control.Sox11[i], 
#              kunche.and.reg$logFC.Control.Sox11MUT[i], 
#              ylim=c(-5, 0), 
#              col=c("red", "green"), notch=TRUE, ylab="log2FC", main = "")
#      dev.off()
#   }
#
#}

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
### we shall make the BARPLOTS instead of BOXPLOTS, 
### and we shall make the plots PER GENE basis ...
##################################################################################################
##################################################################################################
# i = 1

#barplot(c(kunche.and.reg$Control_1.FPKM[i], 
#            kunche.and.reg$Control_2.FPKM[i], 
#            kunche.and.reg$Sox11_1.FPKM[i], 
#            kunche.and.reg$Sox11_2.FPKM[i], 
#            kunche.and.reg$Sox11MUT_1.FPKM[i],
#            kunche.and.reg$Sox11MUT_2.FPKM[i]),
v#            ylim=c(0,50), 
#            main="X",
#            xlab="samples",
#            ylab="gene expression (FPKM)",
#            col=c("blue", "blue", "red", "red", "green", "green") )

# we may use 2^ ..
#barplot(c(kunche.and.reg$logFC.Control.Sox11[i], 
#          kunche.and.reg$logFC.Control.Sox11MUT[i]), 
#          ylim=c(0, 6.5),
#          main="X",
#          xlab="samples",
#          ylab="log2FC", 
#          col=c("red", "green"))

# we may use 2^ ..
#barplot(c(kunche.and.reg$logFC.Control.Sox11[i], 
#          kunche.and.reg$logFC.Control.Sox11MUT[i]), 
#          ylim=c(-5, 0),
#          main="X",
#          xlab="samples",
#          ylab="log2FC", 
#          col=c("red", "green"))

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# PER GENE basis :

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Bmp6"	# DOWN
head(z)

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,20), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-20, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Opn4" #	DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,4), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-20, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Csgalnact1" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,4), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-20, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Spp1" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,40), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-8, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Tgfbr2" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,2), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-8, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Sema3c" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,2), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-8, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Vegfb" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,100), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-8, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Cartpt" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,1600), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-8, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Lsamp" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,40), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-4, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Vegfa" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,30), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-4, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Efna5" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,20), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-4, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Shh" 	# DOWN

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,20), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "DOWN", "png", sep="."), width = 200, height = 600)
barplot(c(-2^(-z$logFC.Control.Sox11), 
          -2^(-z$logFC.Control.Sox11MUT)), 
          ylim=c(-4, 0),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Bcl2" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,10), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Mef2c" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,4), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Dpysl4" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,140), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Plxnb1" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,40), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Ephb3" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,20), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Tfap2b" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,10), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Kif26a" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,4), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,4),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Ngfr" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,70), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,16),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

GENE="Vcan" 	# UP

z <- kunche.and.reg[kunche.and.reg$gene == GENE,] 
head(z)

png(paste(GENE, "the.barplots.of.FPKM", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(z$Control_1.FPKM, 
          z$Control_2.FPKM, 
          z$Sox11_1.FPKM, 
          z$Sox11_2.FPKM, 
          z$Sox11MUT_1.FPKM,
          z$Sox11MUT_2.FPKM),
          ylim=c(0,30), 
          main=GENE,
          xlab="samples",
          ylab="gene expression (FPKM)",
          col=c("blue", "blue", "red", "red", "green", "green") )
dev.off()

# we may use 2^ ..
png(paste(GENE, "the.barplots.of.FC", "UP", "png", sep="."), width = 200, height = 600)
barplot(c(2^(z$logFC.Control.Sox11), 
          2^(z$logFC.Control.Sox11MUT)), 
          ylim=c(0,130),
          main=GENE,
          xlab="samples",
          ylab="FC", 
          col=c("red", "green"))
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###################################################################
################################################################### # to make a VOLCANO PLOT for both UP and DOWN-reg genes :
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# dim(DOWN.moreMUT.FC )
# [1] 173   3
# dim(UP.moreMUT.FC )
# [1] 181   3
# dim(DOWN)
# [1] 902  38
# dim(UP)
# [1] 911  38

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

# marking the DOWN GENES
# DOWN[(DOWN$logFC.Control.Sox11MUT < DOWN$logFC.Control.Sox11),]

DOWN$notes = ""

DOWN$notes[(DOWN$logFC.Control.Sox11MUT < DOWN$logFC.Control.Sox11)] = "down_more_MUT"
# table(DOWN$notes)
#
#              down_more_MUT 
#          729           173 

# dim(DOWN.moreMUT)
# [1] 173  38


########################################################################################
########################################################################################

# marking the UP GENES

UP$notes = ""

UP$notes[(UP$logFC.Control.Sox11MUT > UP$logFC.Control.Sox11)] = "up_more_MUT"
# table(UP$notes)
#
#            up_more_MUT 
#        730         181

# dim(UP.moreMUT)
# [1] 181  38

########################################################################################
########################################################################################
# combining BOTH LISTS :

REG2 <- rbind(DOWN, UP)

dim(REG2)
# [1] 1813   39

write.table(REG2, file="the.list.of.genes.regulated.by.MUT.with.MARKS.txt", sep="\t", quote=F, 
            row.names=F, col.names=T)

# "down_more_MUT"
# "up_more_MUT"

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

library("limma")
library("ggplot2")

#> colnames(REG2)
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
#[39] "notes"

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

png("z.picture.with.genes.regulated.by.MUT.and.COLOR.more.regulated.by.MUTANT.png")
plotWithHighlights(REG2$logFC.Control.Sox11MUT, 
                   -log10(REG2$adj.P.Val.Control.Sox11MUT), 
                   status=REG2$notes,
                   values=c("down_more_MUT", "up_more_MUT"), 
                   bg.col="grey",
                   xlim=c(-8,8), 
                   ylim=c(0,3),
                   hl.cex=0.6, cex.main=0.8, cex.lab =0.8,
                   xlab="log2FC", ylab="-log10 adj.P.Val", 
                   legend= "topright", main="" )
dev.off()

### the display in ggplot2 :

p <- ggplot(REG2, aes(x=logFC.Control.Sox11MUT, 
                 y=-log10(adj.P.Val.Control.Sox11MUT), 
                 col=factor(notes))) +
       geom_point(size=1) +
       theme_bw() +
       xlim(-5, 5) +
       ylim(0, 5) +
       theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
       scale_colour_manual(values = c("grey","down_more_MUT"="green", "up_more_MUT"="red")) +
       labs(x="log2FC", y="-log10 adj.P.Val") + ggtitle("")

png("z.picture.with.genes.regulated.by.MUT.and.COLOR.more.regulated.by.MUTANT.ggplot2.png")
p
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###################################################################################################################################### in order to add the LABELS :
# to integrate REG2 and kunche lists

REG3 <- merge(REG2, kunche, 
              by.x="gene_name", 
              by.y="gene", all.x=TRUE)

dim(REG2)
dim(REG3)

#table(REG3$REGULATION)
#
#DOWN   UP 
#  12    9 

# to prepare the names for the SELECTED GENES

REG3$SELECTED =  ifelse(!is.na(REG3$REGULATION), "SEL", "")
REG3$GENE =  ifelse(!is.na(REG3$REGULATION), REG3$gene_name, "") 

# to write for verification :

write.table(REG3, file="the.list.of.genes.regulated.with.MARKS.and.HIGHLIGHTS.txt", 
            sep="\t", quote=FALSE, row.names=F, col.names=T)


# colnames(REG3)
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
#[39] "notes"                      "REGULATION"                
#[41] "SELECTED"                   "GENE

###################################################################################################################################### 
######################################################################################################################################
##################################################################
##################################################################
### these are the SELECTED GENES :
######################################################################################################################################

q <- ggplot(REG3, aes(x=logFC.Control.Sox11MUT, 
                 y=-log10(adj.P.Val.Control.Sox11MUT), 
                 col=factor(SELECTED))) +
       geom_point(size=1) +
       theme_bw() +
       xlim(-5, 5) +
       ylim(0, 5) +
       theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
       scale_colour_manual(values = c("grey","SEL"="brown")) +
       labs(x="log2FC", y="-log10 adj.P.Val") + ggtitle("")

png("z.picture.with.genes.regulated.by.MUT.and.selected.png", width = 600, height = 600)
q
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
# https://github.com/kevinblighe/EnhancedVolcano
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
# LABELLING THE GENES :
# https://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html
######################################################################################################################################
######################################################################################################################################

q + geom_text(data=REG3, aes(label=GENE))

library(ggrepel)

r = q + geom_text_repel(data=REG3, aes(label=GENE)) +
     ggtitle("selected genes")  

png("z.picture.with.genes.regulated.by.MUT.and.selected.with.LABELS.png")
r
dev.off()

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###################################################################################################################################### to combine the PLOTS :

pqr <- ggplot(REG3, aes(x=logFC.Control.Sox11MUT, 
                 y=-log10(adj.P.Val.Control.Sox11MUT), 
                 col=factor(notes))) +
       geom_point(size=1) +
       theme_bw() +
       xlim(-5, 5) +
       ylim(0, 5) +
       theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
       scale_colour_manual(values = c("grey","down_more_MUT"="green", "up_more_MUT"="red")) +
       labs(x="log2FC", y="-log10 adj.P.Val") +  
       geom_text_repel(data=REG3, aes(label=GENE)) +
       ggtitle("selected genes") 


png("z.picture.with.genes.regulated.by.MUT.and.selected.with.LABELS.integrated.png", width = 800, height = 800)
pqr
dev.off()


pqr
ggsave("z.picture.with.genes.regulated.by.MUT.and.selected.by.with.LABELS.integrated.png")

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###################################################################################################################################### to use ENHANCED VOLCANO

# BiocManager::install('EnhancedVolcano')

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###################################################################################################################################
