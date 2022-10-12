############################################################################################
############################################################################################
############################################################################################
############################################################################################

library("ggplot2")
library("reshape2")
library("data.table")

############################################################################################
############################################################################################
######################################### reading the files with the GENE EXPRESSION COUNTS :

genes <- read.delim("the_GENES.58381_genes.gencode.v28.basic.annotation.28aug2018.txt",
                     sep="\t", header=T, stringsAsFactors=F)

head(genes) 
dim(genes)

genes.dt <- as.data.table(genes)

head(genes.dt) 
dim(genes.dt)

######## to integrate these files : reading the files and changing the names of the columns
############################################################################################
############################################################################################

Aph1 <- read.delim("sample.Aph1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Aph1.simple <- data.frame( Aph1.gene =  Aph1$gene_id,
                           Aph1.count = Aph1$expected_count,
                           Aph1.TPM =   Aph1$TPM,
                           Aph1.FPKM =  Aph1$FPKM, 
                           stringsAsFactors=F)  

head(Aph1) 
dim(Aph1)

head(Aph1.simple)
dim(Aph1.simple)

############################################################################################
############################################################################################

Aph2 <- read.delim("sample.Aph2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Aph2.simple <- data.frame( Aph2.gene =  Aph2$gene_id,
                           Aph2.count = Aph2$expected_count,
                           Aph2.TPM =   Aph2$TPM,
                           Aph2.FPKM =  Aph2$FPKM, 
                           stringsAsFactors=F)  

head(Aph2) 
dim(Aph2)

head(Aph2.simple)
dim(Aph2.simple)

############################################################################################
############################################################################################

Aph3 <- read.delim("sample.Aph3.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Aph3.simple <- data.frame( Aph3.gene =  Aph3$gene_id,
                           Aph3.count = Aph3$expected_count,
                           Aph3.TPM =   Aph3$TPM,
                           Aph3.FPKM =  Aph3$FPKM, 
                           stringsAsFactors=F)  

head(Aph3) 
dim(Aph3)

head(Aph3.simple)
dim(Aph3.simple)

############################################################################################
############################################################################################

Aph_KH7_1 <- read.delim("sample.Aph_KH7_1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Aph_KH7_1.simple <- data.frame( Aph_KH7_1.gene =  Aph_KH7_1$gene_id,
                                Aph_KH7_1.count = Aph_KH7_1$expected_count,
                                Aph_KH7_1.TPM =   Aph_KH7_1$TPM,
                                Aph_KH7_1.FPKM =  Aph_KH7_1$FPKM, 
                                stringsAsFactors=F)  

head(Aph_KH7_1) 
dim(Aph_KH7_1)

head(Aph_KH7_1.simple)
dim(Aph_KH7_1.simple)

############################################################################################
############################################################################################

Aph_KH7_2 <- read.delim("sample.Aph_KH7_2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Aph_KH7_2.simple <- data.frame( Aph_KH7_2.gene =  Aph_KH7_2$gene_id,
                                Aph_KH7_2.count = Aph_KH7_2$expected_count,
                                Aph_KH7_2.TPM =   Aph_KH7_2$TPM,
                                Aph_KH7_2.FPKM =  Aph_KH7_2$FPKM, 
                                stringsAsFactors=F)  

head(Aph_KH7_2) 
dim(Aph_KH7_2)

head(Aph_KH7_2.simple)
dim(Aph_KH7_2.simple)

############################################################################################
############################################################################################

Aph_KH7_3 <- read.delim("sample.Aph_KH7_3.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Aph_KH7_3.simple <- data.frame( Aph_KH7_3.gene =  Aph_KH7_3$gene_id,
                                Aph_KH7_3.count = Aph_KH7_3$expected_count,
                                Aph_KH7_3.TPM =   Aph_KH7_3$TPM,
                                Aph_KH7_3.FPKM =  Aph_KH7_3$FPKM, 
                                stringsAsFactors=F)  

head(Aph_KH7_3) 
dim(Aph_KH7_3)

head(Aph_KH7_3.simple)
dim(Aph_KH7_3.simple)

############################################################################################
############################################################################################

DMSO1_lane1 <- read.delim("sample.DMSO1_lane1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

DMSO1_lane1.simple <- data.frame( DMSO1_lane1.gene =  DMSO1_lane1$gene_id,
                                  DMSO1_lane1.count = DMSO1_lane1$expected_count,
                                  DMSO1_lane1.TPM =   DMSO1_lane1$TPM,
                                  DMSO1_lane1.FPKM =  DMSO1_lane1$FPKM, 
                                  stringsAsFactors=F)  

head(DMSO1_lane1) 
dim(DMSO1_lane1)

head(DMSO1_lane1.simple)
dim(DMSO1_lane1.simple)

############################################################################################
############################################################################################

DMSO1_lane2 <- read.delim("sample.DMSO1_lane2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

DMSO1_lane2.simple <- data.frame( DMSO1_lane2.gene =  DMSO1_lane2$gene_id,
                                  DMSO1_lane2.count = DMSO1_lane2$expected_count,
                                  DMSO1_lane2.TPM =   DMSO1_lane2$TPM,
                                  DMSO1_lane2.FPKM =  DMSO1_lane2$FPKM, 
                                  stringsAsFactors=F)  

head(DMSO1_lane2) 
dim(DMSO1_lane2)

head(DMSO1_lane2.simple)
dim(DMSO1_lane2.simple)

############################################################################################
############################################################################################

DMSO2_lane1 <- read.delim("sample.DMSO2_lane1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

DMSO2_lane1.simple <- data.frame( DMSO2_lane1.gene =  DMSO2_lane1$gene_id,
                                  DMSO2_lane1.count = DMSO2_lane1$expected_count,
                                  DMSO2_lane1.TPM =   DMSO2_lane1$TPM,
                                  DMSO2_lane1.FPKM =  DMSO2_lane1$FPKM, 
                                  stringsAsFactors=F)  

head(DMSO2_lane1) 
dim(DMSO2_lane1)

head(DMSO2_lane1.simple)
dim(DMSO2_lane1.simple)

############################################################################################
############################################################################################

DMSO2_lane2 <- read.delim("sample.DMSO2_lane2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

DMSO2_lane2.simple <- data.frame( DMSO2_lane2.gene =  DMSO2_lane2$gene_id,
                                  DMSO2_lane2.count = DMSO2_lane2$expected_count,
                                  DMSO2_lane2.TPM =   DMSO2_lane2$TPM,
                                  DMSO2_lane2.FPKM =  DMSO2_lane2$FPKM, 
                                  stringsAsFactors=F)  

head(DMSO2_lane2) 
dim(DMSO2_lane2)

head(DMSO2_lane2.simple)
dim(DMSO2_lane2.simple)

############################################################################################
############################################################################################

DMSO3_lane1 <- read.delim("sample.DMSO3_lane1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

DMSO3_lane1.simple <- data.frame( DMSO3_lane1.gene =  DMSO3_lane1$gene_id,
                                  DMSO3_lane1.count = DMSO3_lane1$expected_count,
                                  DMSO3_lane1.TPM =   DMSO3_lane1$TPM,
                                  DMSO3_lane1.FPKM =  DMSO3_lane1$FPKM, 
                                  stringsAsFactors=F)  

head(DMSO3_lane1) 
dim(DMSO3_lane1)

head(DMSO3_lane1.simple)
dim(DMSO3_lane1.simple)

############################################################################################
############################################################################################

DMSO3_lane2 <- read.delim("sample.DMSO3_lane2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

DMSO3_lane2.simple <- data.frame( DMSO3_lane2.gene =  DMSO3_lane2$gene_id,
                                  DMSO3_lane2.count = DMSO3_lane2$expected_count,
                                  DMSO3_lane2.TPM =   DMSO3_lane2$TPM,
                                  DMSO3_lane2.FPKM =  DMSO3_lane2$FPKM, 
                                  stringsAsFactors=F)  

head(DMSO3_lane2) 
dim(DMSO3_lane2)

head(DMSO3_lane2.simple)
dim(DMSO3_lane2.simple)

############################################################################################
############################################################################################

KH7_1 <- read.delim("sample.KH7_1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

KH7_1.simple <- data.frame( KH7_1.gene =  KH7_1$gene_id,
                            KH7_1.count = KH7_1$expected_count,
                            KH7_1.TPM =   KH7_1$TPM,
                            KH7_1.FPKM =  KH7_1$FPKM, 
                                stringsAsFactors=F)  

head(KH7_1) 
dim(KH7_1)

head(KH7_1.simple)
dim(KH7_1.simple)

############################################################################################
############################################################################################

KH7_2 <- read.delim("sample.KH7_2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

KH7_2.simple <- data.frame( KH7_2.gene =  KH7_2$gene_id,
                            KH7_2.count = KH7_2$expected_count,
                            KH7_2.TPM =   KH7_2$TPM,
                            KH7_2.FPKM =  KH7_2$FPKM, 
                                stringsAsFactors=F)  

head(KH7_2) 
dim(KH7_2)

head(KH7_2.simple)
dim(KH7_2.simple)

############################################################################################
############################################################################################

KH7_3 <- read.delim("sample.KH7_3.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

KH7_3.simple <- data.frame( KH7_3.gene =  KH7_3$gene_id,
                            KH7_3.count = KH7_3$expected_count,
                            KH7_3.TPM =   KH7_3$TPM,
                            KH7_3.FPKM =  KH7_3$FPKM, 
                                stringsAsFactors=F)  

head(KH7_3) 
dim(KH7_3)

head(KH7_3.simple)
dim(KH7_3.simple)

############################################################################################
############################################################################################

Noc_1 <- read.delim("sample.Noc_1.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Noc_1.simple <- data.frame( Noc_1.gene =  Noc_1$gene_id,
                            Noc_1.count = Noc_1$expected_count,
                            Noc_1.TPM =   Noc_1$TPM,
                            Noc_1.FPKM =  Noc_1$FPKM, 
                                stringsAsFactors=F)  

head(Noc_1) 
dim(Noc_1)

head(Noc_1.simple)
dim(Noc_1.simple)

############################################################################################
############################################################################################

Noc_2 <- read.delim("sample.Noc_2.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Noc_2.simple <- data.frame( Noc_2.gene =  Noc_2$gene_id,
                            Noc_2.count = Noc_2$expected_count,
                            Noc_2.TPM =   Noc_2$TPM,
                            Noc_2.FPKM =  Noc_2$FPKM, 
                                stringsAsFactors=F)  

head(Noc_2) 
dim(Noc_2)

head(Noc_2.simple)
dim(Noc_2.simple)

############################################################################################
############################################################################################

Noc_3 <- read.delim("sample.Noc_3.rsem.genes.results", sep="\t", header=T, stringsAsFactors=F)

Noc_3.simple <- data.frame( Noc_3.gene =  Noc_3$gene_id,
                            Noc_3.count = Noc_3$expected_count,
                            Noc_3.TPM =   Noc_3$TPM,
                            Noc_3.FPKM =  Noc_3$FPKM, 
                                stringsAsFactors=F)  

head(Noc_3) 
dim(Noc_3)

head(Noc_3.simple)
dim(Noc_3.simple)

############################################################################################
############################################################################################
############################################################################################
############################################################################################

library(data.table)

############################ now integrating these data structures ; we can make DATA TABLES :

Aph1.simple.dt <- as.data.table(Aph1.simple)
Aph2.simple.dt <- as.data.table(Aph2.simple)
Aph3.simple.dt <- as.data.table(Aph3.simple)

Aph_KH7_1.simple.dt <- as.data.table(Aph_KH7_1.simple) 
Aph_KH7_2.simple.dt <- as.data.table(Aph_KH7_2.simple) 
Aph_KH7_3.simple.dt <- as.data.table(Aph_KH7_3.simple) 

DMSO1_lane1.simple.dt <- as.data.table(DMSO1_lane1.simple)
DMSO1_lane2.simple.dt <- as.data.table(DMSO1_lane2.simple)

DMSO2_lane1.simple.dt <- as.data.table(DMSO2_lane1.simple)
DMSO2_lane2.simple.dt <- as.data.table(DMSO2_lane2.simple)

DMSO3_lane1.simple.dt <- as.data.table(DMSO3_lane1.simple)
DMSO3_lane2.simple.dt <- as.data.table(DMSO3_lane2.simple)

KH7_1.simple.dt <- as.data.table(KH7_1.simple)
KH7_2.simple.dt <- as.data.table(KH7_2.simple)
KH7_3.simple.dt <- as.data.table(KH7_3.simple)

Noc_1.simple.dt <- as.data.table(Noc_1.simple)
Noc_2.simple.dt <- as.data.table(Noc_2.simple)
Noc_3.simple.dt <- as.data.table(Noc_3.simple)

############################################################################################
############################################################################################

library(data.table)

setkeyv(genes.dt, c('GENE_ID'))

setkeyv(Aph1.simple.dt, c('Aph1.gene'))
setkeyv(Aph2.simple.dt, c('Aph2.gene'))
setkeyv(Aph3.simple.dt, c('Aph3.gene'))

setkeyv(Aph_KH7_1.simple.dt, c('Aph_KH7_1.gene')) 
setkeyv(Aph_KH7_2.simple.dt, c('Aph_KH7_2.gene'))  
setkeyv(Aph_KH7_3.simple.dt, c('Aph_KH7_3.gene'))  

setkeyv(DMSO1_lane1.simple.dt, c('DMSO1_lane1.gene'))  
setkeyv(DMSO1_lane2.simple.dt, c('DMSO1_lane2.gene'))  

setkeyv(DMSO2_lane1.simple.dt, c('DMSO2_lane1.gene'))  
setkeyv(DMSO2_lane2.simple.dt, c('DMSO2_lane2.gene'))  

setkeyv(DMSO3_lane1.simple.dt, c('DMSO3_lane1.gene'))  
setkeyv(DMSO3_lane2.simple.dt, c('DMSO3_lane2.gene'))  

setkeyv(KH7_1.simple.dt, c('KH7_1.gene'))
setkeyv(KH7_2.simple.dt, c('KH7_2.gene')) 
setkeyv(KH7_3.simple.dt, c('KH7_3.gene'))

setkeyv(Noc_1.simple.dt, c('Noc_1.gene'))
setkeyv(Noc_2.simple.dt, c('Noc_2.gene'))
setkeyv(Noc_3.simple.dt, c('Noc_3.gene'))

############################################################################################
############################################################################################ 
############################################ to integrate ALL the dataframes :

# expression.Aph123 <- genes.dt[Aph1.simple.dt,][Aph2.simple.dt,][Aph3.simple.dt,]

# expression.Aph_KH7_123 <- genes.dt[Aph_KH7_1.simple.dt,][Aph_KH7_2.simple.dt,][Aph_KH7_3.simple.dt,]  

# expression.DMSO <- genes.dt[DMSO1_lane1.simple.dt,][DMSO1_lane2.simple.dt,][DMSO2_lane1.simple.dt,][DMSO2_lane2.simple.dt,][DMSO3_lane1.simple.dt,][DMSO3_lane2.simple.dt,] 

# expression.KH7_123 <- genes.dt[KH7_1.simple.dt,][KH7_2.simple.dt,][KH7_3.simple.dt,] 

# expression.Noc_123 <- genes.dt[Noc_1.simple.dt,][Noc_2.simple.dt,][Noc_3.simple.dt,] 

expression.all.samples <- genes.dt[DMSO1_lane1.simple.dt,][DMSO1_lane2.simple.dt,][DMSO2_lane1.simple.dt,][DMSO2_lane2.simple.dt,][DMSO3_lane1.simple.dt,][DMSO3_lane2.simple.dt,][Aph1.simple.dt,][Aph2.simple.dt,][Aph3.simple.dt,][Aph_KH7_1.simple.dt,][Aph_KH7_2.simple.dt,][Aph_KH7_3.simple.dt,][KH7_1.simple.dt,][KH7_2.simple.dt,][KH7_3.simple.dt,][Noc_1.simple.dt,][Noc_2.simple.dt,][Noc_3.simple.dt,]

expression.all.samples
dim(expression.all.samples)

############################################################################################
############################################################################################
###################### to print the RESULTS, where we have integrated ALL the data frames :

name <- "the_GENES.58381_genes.gencode.v28.basic.annotation.28aug2018.txt"

write.table(expression.all.samples, 
            file=paste(name, ".INTEGRATED.file.ALL.samples.txt", sep=""), 
            sep="\t", quote=FALSE, 
            row.names = FALSE, col.names = TRUE)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
