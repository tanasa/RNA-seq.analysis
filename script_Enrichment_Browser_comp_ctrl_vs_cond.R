######################################################################
######################################################################
######################################################################

library("EnrichmentBrowser")
library("SummarizedExperiment")
library("ReportingTools")

library(stringr)
library(dplyr)
library(tidyr)

######################################################################
######################################################################
######################################################################
######################################################################

NAME = "compare.CTRL.COND"

file.genes = read.delim("analysis.LIMMA.integrating.all.samples.with.data.table.printing.FPKM.and.CELL.CYCLE.txt", 
                         header=T, sep="\t", stringsAsFactors=F)


file.genes.subset = subset(file.genes, select = c("GENE_ID", 
                                                  "GENE_NAME", "GENE_NAME_ID",

"DMSO1_lane1.count",      
"DMSO2_lane1.count",
"DMSO3_lane1.count",  

"DMSO1_lane2.count",
"DMSO2_lane2.count",
"DMSO3_lane2.count",            

"Aph1.count",             
"Aph2.count",  
"Aph3.count",   

"Aph_KH7_1.count",       
"Aph_KH7_2.count",       
"Aph_KH7_3.count", 

"KH7_1.count",            
"KH7_2.count", 
"KH7_3.count",            

"Noc_1.count",                          
"Noc_2.count",            
"Noc_3.count",

"logFC.Aph",              
"P.Value.Aph",
"adj.P.Val.Aph",          
               
"logFC.Aph_KH7",
"P.Value.Aph_KH7",         
"adj.P.Val.Aph_KH7",      
          
"logFC.KH7",              
"P.Value.KH7",
"adj.P.Val.KH7",          
             
"logFC.Noc",              
"P.Value.Noc",
"adj.P.Val.Noc" )) 

######################################################################
######################################################################
######################################################################

write.table(file.genes.subset, 
            file = paste(NAME, "the.subsets.txt", sep="."), 
            sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE) 

######################################################################
######################################################################
######################################################################

file.genes.EXP = subset(file.genes, select = c("GENE_ID", 
                                                  
"DMSO1_lane1.count",      
"DMSO2_lane1.count",
"DMSO3_lane1.count",  

"Aph1.count",       
"Aph2.count",       
"Aph3.count", 
       
"logFC.Aph",
"P.Value.Aph",         
"adj.P.Val.Aph")) 

######################################################################
######################################################################
###################################################################### 

library(stringr)
library(dplyr)
library(tidyr)

file.genes.EXP = file.genes.EXP %>% tidyr::separate(GENE_ID, c("ENSEMBL"), sep="\\.", extra="drop")

file.genes.EXP = file.genes.EXP %>% distinct()

######################################################################
######################################################################
######################################################################

dim(file.genes.EXP)      # [1] 12944    10
colnames(file.genes.EXP)

######################################################################
######################################################################
###################################################################### 

file.genes.EXP$FC         = file.genes.EXP$logFC.Aph  
file.genes.EXP$PVAL       = file.genes.EXP$P.Value.Aph       
file.genes.EXP$ADJ.PVAL   = file.genes.EXP$adj.P.Val.Aph

######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### to set up the COUNTS

counts = subset(file.genes.EXP, select = c("ENSEMBL", 
                                            "DMSO1_lane1.count",      
                                            "DMSO2_lane1.count",
                                            "DMSO3_lane1.count",  

                                            "Aph1.count",       
                                            "Aph2.count",       
                                            "Aph3.count" ))

rownames(counts) = counts$ENSEMBL

counts = counts[,-1]

counts = as.matrix(counts)

head(counts)

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### to set up the colData

rowData = subset(file.genes.EXP, select = c("ENSEMBL", "FC", "PVAL", "ADJ.PVAL")) 

rownames(rowData) = rowData$ENSEMBL

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### to set up the rowData

colData = data.frame( SAMPLE = c("DMSO1_lane1.count", "DMSO2_lane1.count", "DMSO3_lane1.count",  
                                 "Aph1.count", "Aph2.count", "Aph3.count" )) 

colData$GROUP = c(0,0,0,1,1,1)

###################################################################### RNA-seq :
######################################################################
######################################################################

DATA = SummarizedExperiment(assays=list(counts=counts),
                             rowData=rowData, 
                             colData=colData)

###################################################################### RangedSummarizedExperiment
######################################################################
###################################################################### 

rowData(DATA)
colData(DATA)

(DATA)

rownames(DATA)
colnames(DATA)

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

png(paste(NAME, "the.p.values.distributions.and.a.volcano.plot", "png", sep="."))
par(mfrow = c(1,2))
pdistr(rowData(DATA)$PVAL)
volcano(rowData(DATA)$FC, rowData(DATA)$ADJ.PVAL)
dev.off()

#############################################################################
#############################################################################
############################################################################# computing the differential expression
# only for curiosity to recompute the DIFFERENTIAL EXPRESSION

DATA2 = DATA

rowData(DATA2, use.names=TRUE)
colData(DATA2, use.names=TRUE)

#############################################################################
#############################################################################
#############################################################################
############################################################################# 

idTypes("hsa")

# idTypes("hsa")
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
#[11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
#[16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
#[21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
#[26] "UNIPROT"

# idTypes("mmu")
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MGI"          "ONTOLOGY"    
# [16] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"     
# [21] "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"

# ID mapping for the airway dataset (from ENSEMBL to ENTREZ gene ids) can then be carried out using the function idMap.

DATA <- idMap(DATA, org = "hsa", from = "ENSEMBL", to = "ENTREZID")

#############################################################################
#############################################################################

head(colData(DATA))
head(assays(DATA)$counts)
head(rowData(DATA))

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

# db: Database from which gene sets should be retrieved. Currently,
#             either 'go' (default), 'kegg', 'msigdb', or 'enrichr'.

#   cache: Logical.  Should a locally cached version used if available?
#          Defaults to 'TRUE'.

# return.type: Character. Determines whether gene sets are returned as a
#          simple list of gene sets (each being a character vector of
#          gene IDs), or as an object of class 'GeneSetCollection'.

#     ...: Additional arguments for individual gene set databases.  For
#          'db = "GO"':

#            • onto: Character. Specifies one of the three GO
#              ontologies: 'BP' (biological process), 'MF' (molecular
#              function), 'CC' (cellular component). Defaults to 'BP'.

#            • mode: Character. Determines in which way the gene sets
#              are retrieved. This can be either 'GO.db' or 'biomart'.
#              The 'GO.db' mode creates the gene sets based on BioC
#              annotation packages - which is fast, but represents not
#              necessarily the most up-to-date mapping. In addition,
#              this option is only available for the currently supported
#              model organisms in BioC.  The 'biomart' mode downloads
#              the mapping from BioMart - which can be time consuming,
#              but allows to select from a larger range of organisms and
#             contains the latest mappings.  Defaults to 'GO.db'.

#          For 'db = "msigdb":'
#
#            • cat: Character.  MSigDB collection category: 'H'
#              (hallmark), 'C1' (genomic position), 'C2' (curated
#              databases), 'C3' (binding site motifs), 'C4'
#              (computational cancer), 'C5' (Gene Ontology), 'C6'
#              (oncogenic), 'C7' (immunologic). See references.

#            • subcat: Character. MSigDB collection subcategory. Depends
#              on the chosen MSigDB collection category. For example,
#              'MIR' to obtain microRNA targets from the 'C3'
#              collection. See references.

#          For 'db = "enrichr"':

#            • lib: Character. Enrichr gene set library. For example,
#              'Genes_Associated_with_NIH_Grants' to obtain gene sets
#              based on associations with NIH grants. See references.

# Examples:
#
#         # (1) Typical usage for gene set enrichment analysis with GO:
#         # Biological process terms based on BioC annotation (for human)
#         go.gs <- getGenesets(org="hsa", db="go")
         
         # eq.:  
         # go.gs <- getGenesets(org="hsa", db="go", onto="BP", mode="GO.db")
         
         # Alternatively:
         # downloading from BioMart 
         # this may take a few minutes ...
#         go.gs <- getGenesets(org="hsa", db="go", mode="biomart")
         
         # (2) Defining gene sets according to KEGG  
#         kegg.gs <- getGenesets(org="hsa", db="kegg")
         
         # (3) Obtaining *H*allmark gene sets from MSigDB
#         hall.gs <- getGenesets(org="hsa", db="msigdb", cat="H")
         
         
         # (4) Obtaining gene sets from Enrichr
#         tfppi.gs <- getGenesets(org="hsa", db="enrichr", lib="Transcription_Factor_PPIs")
         
         # displaying available Enrichr gene set libraries
#         getGenesets(org="hsa", db="enrichr", show.libs=TRUE)        
         
         # (6) parsing gene sets from GMT
#         gmt.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
#         gs <- getGenesets(gmt.file)     
         
         # (7) writing gene sets to file
#         writeGMT(gs, gmt.file)

#############################################################################
#############################################################################
#############################################################################

kegg.gs <- getGenesets(org="hsa", db="kegg")

go.gs <- getGenesets(org="hsa", db="go", go.onto="BP", go.mode="GO.db")

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# SBEA

# 'ora': overrepresentation analysis, simple and frequently used
#     test based on the hypergeometric distribution (see Goeman and
#     Buhlmann, 2007, for a critical review).

#     'safe': significance analysis of function and expression,
#     generalization of ORA, includes other test statistics, e.g.
#     Wilcoxon's rank sum, and allows to estimate the significance of
#     gene sets by sample permutation; implemented in the safe package
#     (Barry et al., 2005).

#     'gsea': gene set enrichment analysis, frequently used and widely
#     accepted, uses a Kolmogorov-Smirnov statistic to test whether the
#     ranks of the p-values of genes in a gene set resemble a uniform
#     distribution (Subramanian et al., 2005).

#     'padog': pathway analysis with down-weighting of overlapping
#     genes, incorporates gene weights to favor genes appearing in few
#     pathways versus genes that appear in many pathways; implemented in
#    the PADOG package.

#     'roast': rotation gene set test, uses rotation instead of
#     permutation for assessment of gene set significance; implemented
#     in the limma and edgeR packages for microarray and RNA-seq data,
#     respectively.

#     'camera': correlation adjusted mean rank gene set test, accounts
#     for inter-gene correlations as implemented in the limma and edgeR
#    packages for microarray and RNA-seq data, respectively.

#     'gsa': gene set analysis, differs from GSEA by using the maxmean
#     statistic, i.e. the mean of the positive or negative part of gene
#     scores in the gene set; implemented in the GSA package.

#     'gsva': gene set variation analysis, transforms the data from a
#     gene by sample matrix to a gene set by sample matrix, thereby
#     allowing the evaluation of gene set enrichment for each sample;
#     implemented in the GSVA package.

#     'globaltest': global testing of groups of genes, general test of
#     groups of genes for association with a response variable;
#     implemented in the globaltest package.

#     'samgs': significance analysis of microarrays on gene sets,
#     extends the SAM method for single genes to gene set analysis (Dinu
#     et al., 2007).

#     'ebm': empirical Brown's method, combines $p$-values of genes in a
#     gene set using Brown's method to combine $p$-values from dependent
#     tests; implemented in the EmpiricalBrownsMethod package.

#     'mgsa': model-based gene set analysis, Bayesian modeling approach
#     taking set overlap into account by working on all sets
#     simultaneously, thereby reducing the number of redundant sets;
#    implemented in the mgsa package.

#     It is also possible to use additional set-based enrichment
#     methods.  This requires to implement a function that takes 'se',
#     'gs', 'alpha', and 'perm' as arguments and returns a numeric
#     vector 'ps' storing the resulting p-value for each gene set in
#     'gs'. This vector must be named accordingly (i.e. names(ps) ==
#     names(gs)). See examples.

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# NBEA

# nbea(method = EnrichmentBrowser::nbeaMethods(), se, gs, grn,
#       prune.grn = TRUE, alpha = 0.05, perm = 1000,
#       padj.method = "none", out.file = NULL, browse = FALSE, ...)

#  'ggea': gene graph enrichment analysis, scores gene sets according
#     to consistency within the given gene regulatory network, i.e.
#     checks activating regulations for positive correlation and
#     repressing regulations for negative correlation of regulator and
#     target gene expression (Geistlinger et al., 2011). When using
#     'ggea' it is possible to estimate the statistical significance of
#     the consistency score of each gene set in two different ways: (1)
#     based on sample permutation as described in the original
#     publication (Geistlinger et al., 2011) or (2) using an
#     approximation in the spirit of Bioconductor's npGSEA package that
#     is much faster.

#     'spia': signaling pathway impact analysis, combines ORA with the
#     probability that expression changes are propagated across the
#     pathway topology; implemented in Bioconductor's SPIA package
#     (Tarca et al., 2009).

#     'pathnet': pathway analysis using network information, applies ORA
#     on combined evidence for the observed signal for gene nodes and
#     the signal implied by connected neighbors in the network;
#     implemented in Bioconductor's PathNet package.

#    'degraph': differential expression testing for gene graphs,
#     multivariate testing of differences in mean incorporating
#     underlying graph structure; implemented in Bioconductor's DEGraph
#     package.

#     'topologygsa': topology-based gene set analysis, uses Gaussian
#     graphical models to incorporate the dependence structure among
#     genes as implied by pathway topology; implemented in CRAN's
#     topologyGSA package.

#     'ganpa': gene association network-based pathway analysis,
#     incorporates network-derived gene weights in the enrichment
#     analysis; implemented in CRAN's GANPA package.

#     'cepa': centrality-based pathway enrichment, incorporates network
#     centralities as node weights mapped from differentially expressed
#     genes in pathways; implemented in CRAN's CePa package.

#     'netgsa': network-based gene set analysis, incorporates external
#     information about interactions among genes as well as novel
#     interactions learned from data; implemented in CRAN's NetGSA
#     package.

#     It is also possible to use additional network-based enrichment
#     methods. This requires to implement a function that takes 'se',
#     'gs', 'grn', 'alpha', and 'perm' as arguments and returns a
#     numeric matrix 'res.tbl' with a mandatory column named 'PVAL'
#     storing the resulting p-value for each gene set in 'gs'. The rows
#     of this matrix must be named accordingly (i.e. rownames(res.tbl)
#     == names(gs)). See examples.

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# ORA

ora <- sbea(method="ora", se=DATA, gs=kegg.gs, 
                perm = 0, 
                # alpha = 0.05,  
                # padj.method = "none", 
                out.file = paste(NAME, "analysis.KEGG.ORA.txt", sep="."))

gsRanking(ora)

write.table(as.data.frame(gsRanking(ora)),
            file = paste(NAME, "analysis.KEGG.ORA.v2.txt", sep="."), 
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(ora)) ) { 
eaBrowse(ora, out.dir= paste(getwd(), "ORA.KEGG", sep="/"), report.name = "analysis.ORA.KEGG")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# GSEA

gsea <- sbea(method="gsea", se=DATA, gs=kegg.gs, 
             # perm = 1000, 
             alpha = 0.2,   ### changing the DEFAULT SETTINGS for GSEA
             # padj.method = "none",
             out.file = paste(NAME, "analysis.KEGG.GSEA.txt", sep="."))  

gsRanking(gsea)

write.table(as.data.frame(gsRanking(gsea)),
            file = paste(NAME, "analysis.KEGG.GSEA.v2.txt", sep="."),
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(gsea)) ) { 
eaBrowse(gsea, out.dir= paste(getwd(), "GSEA.KEGG", sep="/"), report.name = "analysis.GSEA.KEGG")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# ROAST

roast <- sbea(method="roast", se=DATA, gs=kegg.gs,
              # perm=1000, 
              # alpha = 0.05, 
              # padj.method = "none", 
              out.file = paste(NAME, "analysis.KEGG.ROAST.txt", sep="."))

gsRanking(roast) 
 
write.table(as.data.frame(gsRanking(roast)),
            file = paste(NAME, "analysis.KEGG.ROAST.v2.txt", sep="."),
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(roast)) ) { 
eaBrowse(roast,  out.dir = paste(getwd(), "ROAST.KEGG", sep="/"), report.name = "analysis.ROAST.KEGG")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# COMPUTING the NETWORKS :

kegg.grn <- compileGRN(org="hsa", db="kegg")

# downloadPathways(org, cache = TRUE, out.dir = NULL, zip = FALSE)

# compileGRN(org, db = "kegg", act.inh = TRUE, map2entrez = TRUE,
#       keep.type = FALSE, kegg.native = FALSE)
     
# Arguments:
#
#     org: An organism in KEGG three letter code, e.g. 'hsa' for 'Homo
#          sapiens'.  Alternatively, and mainly for backward
#          compatibility, this can also be either a list of
#          'KEGGPathway' objects or an absolute file path of a zip
#          compressed archive of pathway xml files in KGML format.

#      db: Pathway database.  This should be one or more DBs out of
#          'kegg', 'reactome', 'biocarta', and 'nci'.  See
#          'pathwayDatabases' for available DBs of the respective
#           organism.  Default is 'kegg'

#     'pathwayDatabases', 'pathways', 'KEGGPathway', 'parseKGML',
#     'downloadPathways'

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# SPIA

spia <- nbea(method="spia", se=DATA, gs=kegg.gs, grn=kegg.grn, 
             # perm=1000, 
             alpha = 0.2,               ###### 've changed the p-value to 0.2 !!!! 
             # padj.method = "none",
              out.file = paste(NAME, "analysis.KEGG.SPIA.txt", sep="."))

gsRanking(spia)

write.table(as.data.frame(gsRanking(spia)),
            file = paste(NAME, "analysis.KEGG.SPIA.v2.txt", sep="."), 
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(spia)) ) { 
eaBrowse(spia,  out.dir = paste(getwd(), "SPIA.KEGG", sep="/"), report.name = "analysis.SPIA.KEGG")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# GGEA

ggea <- nbea(method="ggea", se=DATA, gs=kegg.gs, grn=kegg.grn, 
             # perm=1000, 
             alpha = 0.2,           ###### I've changed the p-value to 0.2 ! 
             # padj.method = "none",
             out.file = paste(NAME, "analysis.KEGG.GGEA.txt", sep="."))

gsRanking(ggea)

write.table(as.data.frame(gsRanking(ggea)),
            file = paste(NAME, "analysis.KEGG.GGEA.v2.txt", sep="."),
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(ggea)) ) { 
eaBrowse(ggea,  out.dir = paste(getwd(), "GGEA.KEGG", sep="/"), report.name = "analysis.GGEA.KEGG")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# GGEA GRAPH

png("the.GGEA,graph.EXAMPLE.cell.cycle.png",  width = 100, height = 100, units = "cm", res=100)
par(mfrow=c(1,2))

ggeaGraph(gs=kegg.gs[["hsa04110_Cell_cycle"]],
          grn=kegg.grn, 
          se=DATA)

ggeaGraphLegend()

dev.off()

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# ORA

ora2 <- sbea(method="ora", se=DATA, gs=go.gs, 
                perm=0, 
                # alpha = 0.05,  
                # padj.method = "none", 
                out.file = paste(NAME, "analysis.GO.ORA.txt", sep="."))

gsRanking(ora2)

write.table(as.data.frame(gsRanking(ora2)),
            file = paste(NAME, "analysis.GO.ORA.v2.txt", sep="."), 
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(ora2)) ) { 
eaBrowse(ora2, out.dir= paste(getwd(), "ORA.GO", sep="/"), report.name = "analysis.ORA.GO")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################

# gsea2 <- sbea(method="gsea", se=DATA, gs=go.gs, 
#             # perm=1000, 
#             alpha = 0.2,   ### changing the DEFAULT SETTINGS for GSEA!!!!!! 
#             # padj.method = "none",
#             out.file = paste(NAME, "analysis.GO.GSEA.txt", sep="."))  

# gsRanking(gsea2)

# write.table(as.data.frame(gsRanking(gsea2)),
#            file = paste(NAME, "analysis.GO.GSEA.v2.txt", sep="."),
#            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

# if (!is.null(gsRanking(gsea2)) ) { 
# eaBrowse(gsea2, out.dir= paste(getwd(), "GO.GSEA", sep="/"), report.name = "analysis.GO.GSEA")
# }

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# ROAST

roast2 <- sbea(method="roast", se=DATA, gs=go.gs,
              # perm=1000, 
              # alpha = 0.05, 
              # padj.method = "none", 
              out.file = paste(NAME, "analysis.GO.ROAST.txt", sep="."))

gsRanking(roast2) 
 
write.table(as.data.frame(gsRanking(roast2)),
            file = paste(NAME, "analysis.GO.ROAST.v2.txt", sep="."),
            sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE) 

if (!is.null(gsRanking(roast2)) ) { 
eaBrowse(roast2,  out.dir = paste(getwd(), "ROAST.GO", sep="/"), report.name = "analysis.ROAST.GO")
}

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# COMBINING the LISTS : 

# res.list <- list(ora, gsea)
# comb.res <- combResults(res.list)

# gsRanking(comb.res)

# eaBrowse(comb.res, graph.view=hsa.grn, nr.show=5)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# 

#############################################################################
#############################################################################
#############################################################################
#############################################################################
############################################################################# MSIGDB now version 7.1
# https://www.gsea-msigdb.org/gsea/msigdb

# getGenesets         package:EnrichmentBrowser          R Documentation
# Definition of gene sets according to different sources
# Functionality for retrieving gene sets for an organism under
#     investigation from databases such as GO and KEGG. Parsing and
#     writing a list of gene sets from/to a flat text file in GMT format
#     is also supported.

#     The GMT (Gene Matrix Transposed) file format is a tab delimited
#     file format that describes gene sets.  In the GMT format, each row
#     represents a gene set.  Each gene set is described by a name, a
#     description, and the genes in the gene set. See references.

# getGenesets(org, db = c("go", "kegg"), cache = TRUE,
#       go.onto = c("BP", "MF", "CC"), go.mode = c("GO.db", "biomart"),
#       return.type = c("list", "GeneSetCollection"))
     
#     writeGMT(gs, gmt.file)


# Note #1: Use getGenesets with db = "msigdb" to obtain gene set collections for 11 different
# species from the Molecular Signatures Database (MSigDB). Analogously, getGenesets with
# db = "enrichr" allows to obtain gene set libraries from the comprehensive Enrichr collection
# for 5 different species.

# Note #2: The idMap function can be used to map gene sets from NCBI Entrez Gene IDs to
# other common gene ID types such as ENSEMBL gene IDs or HGNC symbols as described in Section 6.

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### SBEA()
######################################################################

# sbeaMethods()
     
#     sbea(method = EnrichmentBrowser::sbeaMethods(), se, gs, alpha = 0.05,
#       perm = 1000, padj.method = "none", out.file = NULL,
#       browse = FALSE, ...)
     
# method: Set-based enrichment analysis method.  Currently, the
#          following set-based enrichment analysis methods are
#          supported: ‘ora’, ‘safe’, ‘gsea’, ‘padog’, ‘roast’, ‘camera’,
#          ‘gsa’, ‘gsva’, ‘globaltest’, ‘samgs’, ‘ebm’, and ‘mgsa’.  For
#          basic ora also set 'perm=0'. Default is ‘ora’.  This can also
#          be the name of a user-defined function implementing set-based
#          enrichment. See Details.

# se: Expression dataset.  An object of class
#          ‘SummarizedExperiment’.  Mandatory minimal annotations:

#            • colData column storing binary group assignment (named "GROUP")

#            • rowData column storing (log2) fold changes of
#              differential expression between sample groups (named "FC")

#            • rowData column storing adjusted (corrected for multiple
#              testing) p-values of differential expression between
#              sample groups (named "ADJ.PVAL")

#          Additional optional annotations:

#            • colData column defining paired samples or sample blocks
#              (named "BLOCK")

#            • metadata slot named "annotation" giving the organism
#              under investigation in KEGG three letter code (e.g. "hsa"
#              for Homo sapiens)

#            • metadata slot named "dataType" indicating the expression
#              data type ("ma" for microarray, "rseq" for RNA-seq)

# gs: Gene sets.  Either a list of gene sets (character vectors of
#          gene IDs) or a text file in GMT format storing all gene sets
#          under investigation.

# perm: Number of permutations of the sample group assignments.
#          Defaults to 1000. For basic ora set 'perm=0'.  Using
#          method="gsea" and 'perm=0' invokes the permutation
#          approximation from the npGSEA package.

# padj.method: Method for adjusting nominal gene set p-values to multiple
#          testing.  For available methods see the man page of the stats
#          function ‘p.adjust’.  Defaults to'none', i.e. leaves the
#          nominal gene set p-values unadjusted.

# out.file: Optional output file the gene set ranking will be written to.

#    'ora': overrepresentation analysis, simple and frequently used
#     test based on the hypergeometric distribution (see Goeman and
#     Buhlmann, 2007, for a critical review).

#     'safe': significance analysis of function and expression,
#     generalization of ORA, includes other test statistics, e.g.
#     Wilcoxon's rank sum, and allows to estimate the significance of
#     gene sets by sample permutation; implemented in the safe package
#     (Barry et al., 2005).

#     'gsea': gene set enrichment analysis, frequently used and widely
#     accepted, uses a Kolmogorov-Smirnov statistic to test whether the
#     ranks of the p-values of genes in a gene set resemble a uniform
#     distribution (Subramanian et al., 2005).

#     'padog': pathway analysis with down-weighting of overlapping
#     genes, incorporates gene weights to favor genes appearing in few
#     pathways versus genes that appear in many pathways; implemented in
#     the PADOG package.

#     'roast': rotation gene set test, uses rotation instead of
#     permutation for assessment of gene set significance; implemented
#     in the limma and edgeR packages for microarray and RNA-seq data,
#     respectively.

#    'camera': correlation adjusted mean rank gene set test, accounts
#     for inter-gene correlations as implemented in the limma and edgeR
#     packages for microarray and RNA-seq data, respectively.

#     'gsa': gene set analysis, differs from GSEA by using the maxmean
#     statistic, i.e. the mean of the positive or negative part of gene
#     scores in the gene set; implemented in the GSA package.

#     'gsva': gene set variation analysis, transforms the data from a
#     gene by sample matrix to a gene set by sample matrix, thereby
#     allowing the evaluation of gene set enrichment for each sample;
#     implemented in the GSVA package.

#     'globaltest': global testing of groups of genes, general test of
#     groups of genes for association with a response variable;
#     implemented in the globaltest package.

#     'samgs': significance analysis of microarrays on gene sets,
#     extends the SAM method for single genes to gene set analysis (Dinu
#     et al., 2007)

#     'ebm': empirical Brown's method, combines $p$-values of genes in a
#     gene set using Brown's method to combine $p$-values from dependent
#     tests; implemented in the EmpiricalBrownsMethod package.

#     'mgsa': model-based gene set analysis, Bayesian modeling approach
#     taking set overlap into account by working on all sets
#     simultaneously, thereby reducing the number of redundant sets;
#     implemented in the mgsa package.

#     It is also possible to use additional set-based enrichment
#     methods.  This requires to implement a function that takes 'se',
#     'gs', 'alpha', and 'perm' as arguments and returns a numeric
#     vector 'ps' storing the resulting p-value for each gene set in
#     'gs'. This vector must be named accordingly (i.e. names(ps) ==
#     names(gs)).

######################################################################
######################################################################
###################################################################### GO/KEGG overrepresentation analysis
######################################################################
######################################################################
######################################################################
###################################################################### idMap

#  Functionality to map the rownames of a SummarizedExperiment
#     between common gene ID types such as ENSEMBL and ENTREZ.

#  idMap(se, org = NA, from = "ENSEMBL", to = "ENTREZID",
#       multi.to = "first", multi.from = "first")
     
#     idTypes(org)
     
# se: An object of class 'SummarizedExperiment'. Expects the names
#          to be of gene ID type given in argument 'from'.

#     org: Character. Organism in KEGG three letter code, e.g. 'hsa' for
#          'Homo sapiens'.  See references.

#    from: Character. Gene ID type from which should be mapped.
#          Corresponds to the gene ID type of the names of argument
#          'se'. Note that 'from' is ignored if 'to' is a 'rowData'
#          column of 'se'.  Defaults to 'ENSEMBL'.

#      to: Character. Gene ID type to which should be mapped.
#          Corresponds to the gene ID type the rownames of argument 'se'
#          should be updated with. Note that this can also be the name
#         of a column in the 'rowData' slot of 'se' to specify
#         of a column in the 'rowData' slot of 'se' to specify
#          user-defined mappings in which conflicts have been manually
#          resolved. Defaults to 'ENTREZID'.

# multi.to: How to resolve 1:many mappings, i.e. multiple to.IDs for a
#          single from.ID? This is passed on to the 'multiVals' argument
#          of 'mapIds' and can thus take several pre-defined values, but
#          also the form of a user-defined function. However, note that
#          this requires that a single to.ID is returned for each
#          from.ID. Default is '"first"', which accordingly returns the
#          first to.ID mapped onto the respective from.ID.

######################################################################
######################################################################
######################################################################
###################################################################### 
#####################################################################

# sbeaMethods()
#  [1] "ora"        "safe"       "gsea"       "gsa"        "padog"    
#  [6] "globaltest" "roast"      "camera"     "gsva"       "samgs"    
# [11] "ebm"        "mgsa"

# nbeaMethods()
# [1] "ggea"        "spia"        "pathnet"     "degraph"     "ganpa"      
# [6] "cepa"        "topologygsa" "netgsa" 

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### Network-based enrichment analysis

# Compilation of a gene regulatory network from pathway databases
# Description:

#     To perform network-based enrichment analysis a gene regulatory
#     network (GRN) is required. There are well-studied processes and
#     organisms for which comprehensive and well-annotated regulatory
#     networks are available, e.g. the RegulonDB for E. coli and
#     Yeastract for S. cerevisiae.  However, in many cases such a
#     network is missing.  A first simple workaround is to compile a
#     network from regulations in pathway databases such as KEGG.

# Usage:
# compileGRN(org, db = "kegg", act.inh = TRUE, map2entrez = TRUE,
#       keep.type = FALSE, kegg.native = FALSE)

# org: An organism in KEGG three letter code, e.g. ‘hsa’ for ‘Homo
#          sapiens’.  Alternatively, and mainly for backward
#          compatibility, this can also be either a list of
#          ‘KEGGPathway’ objects or an absolute file path of a zip
#          compressed archive of pathway xml files in KGML format.

# db: Pathway database.  This should be one or more DBs out of
#          'kegg', 'reactome', 'biocarta', and 'nci'.  See
#          ‘pathwayDatabases’ for available DBs of the respective
#          organism.  Default is 'kegg'. Note: when dealing with
#          non-model organisms, GRN compilation is currently only
#          supported directly from KEGG (the argument ‘kegg.native’
#          should accordingly be set to ‘TRUE’).

# act.inh: Should gene regulatory interactions be classified as
#          activating (+) or inhibiting (-)?  If TRUE, this will drop
#          interactions for which such a classification cannot be made
#          (e.g. binding events). Otherwise, all interactions found in
#          the pathway DB will be included. Default is ‘TRUE’.

# map2entrez: Should gene identifiers be mapped to NCBI Entrez Gene IDs?
#           This only applies to Reactome and NCI as they both use
#           UNIPROT IDs.  This is typically recommended when using the
#          GRN for network-based enrichment analysis with the
#          EnrichmentBrowser.  Default is ‘TRUE’.

# Having found gene sets that show enrichment for differential expression,
# we are now interested whether these findings can be supported by known regulatory interactions.

#############################################################################
############################################################################# eaBrowse()

# Functions to extract a flat gene set ranking from an enrichment
#     analysis result object and to detailedly explore it.

# eaBrowse(res, nr.show = -1, graph.view = NULL, html.only = FALSE,
#       out.dir = NULL, report.name = NULL)
     
#     gsRanking(res, signif.only = TRUE)
     
# res: Enrichment analysis result list (as returned by the functions
#          'sbea' and 'nbea').

# nr.show: Number of gene sets to show.  As default all statistically
#          significant gene sets are displayed.

# graph.view: Optional.  Should a graph-based summary (reports and
#          visualizes consistency of regulations) be created for the
#          result?  If specified, it needs to be a gene regulatory
#          network, i.e. either an absolute file path to a tabular file
#          or a character matrix with exactly *THREE* cols; 1st col =
#          IDs of regulating genes; 2nd col = corresponding regulated
#          genes; 3rd col = regulation effect; Use '+' and '-' for
#          activation/inhibition.

#############################################################################
############################################################################# nbeaMethods()

# nbeaMethods         package:EnrichmentBrowser          R Documentation
# Network-based enrichment analysis (NBEA)

# This is the main function for network-based enrichment analysis.
#     It implements and wraps existing implementations of several
#     frequently used methods and allows a flexible inspection of
#     resulting gene set rankings.

# nbeaMethods()
     
#     nbea(method = EnrichmentBrowser::nbeaMethods(), se, gs, grn,
#       prune.grn = TRUE, alpha = 0.05, perm = 1000,
#       padj.method = "none", out.file = NULL, browse = FALSE, ...)
     
# method: Network-based enrichment analysis method.  Currently, the
#          following network-based enrichment analysis methods are
#          supported: 'ggea', 'spia', 'pathnet', 'degraph',
#          'topologygsa', 'ganpa', 'cepa', and 'netgsa'. Default is
#          'ggea'.  This can also be the name of a user-defined function
#          implementing network-based enrichment. See Details.

#      se: Expression dataset.  An object of class
#          'SummarizedExperiment'.  Mandatory minimal annotations:

#  <pre> • colData column storing binary group assignment (named GROUP)

# • rowData column storing (log2) fold changes of
#              differential expression between sample groups (named FC)

#            • rowData column storing adjusted (corrected for multiple
#              testing) p-values of differential expression between
#            sample groups (named &quot;ADJ.PVAL&quot;)

#          Additional optional annotations:

#            • colData column defining paired samples or sample blocks
#              (named &quot;BLOCK&quot;)

#            • metadata slot named &quot;annotation&quot; giving the organism
#              under investigation in KEGG three letter code (e.g. &quot;hsa&quot;
#              for Homo sapiens)

#            • metadata slot named &quot;dataType&quot; indicating the expression
#              data type (&quot;ma&quot; for microarray, &quot;rseq&quot; for RNA-seq)

######################################################################
######################################################################
######################################################################
######################################################################
###################################################################### COMPARISONS :

# The function combResults implements the straightforward combination of results,
# thereby facilitating seamless comparison of results across methods.
# For demonstration, we use the ORA and GSEA results for the ALL dataset
# from the previous section:

# combResults         package:EnrichmentBrowser          R Documentation
# Combining enrichment analysis results

# Different enrichment analysis methods usually result in different
#     gene set rankings for the same dataset.  This function allows to
#     combine results from the different set-based and network-based
#     enrichment analysis methods.  This includes the computation of
#     average gene set ranks across methods.

# combResults(res.list, rank.col = configEBrowser("PVAL.COL"),
#       decreasing = FALSE, rank.fun = c("comp.ranks", "rel.ranks",
#       "abs.ranks"), comb.fun = c("mean", "median", "min", "max", "sum"))
     
# res.list: A list of enrichment analysis result lists (as returned by
#          the functions 'sbea' and 'nbea').

# rank.col: Rank column.  Column name of the enrichment analysis result
#          table that should be used to rank the gene sets.  Defaults to
#          the gene set p-value column, i.e. gene sets are ranked
#          according to gene set significance.

# res.list <- list(ora.all, gsea.all)
# comb.res <- combResults(res.list)
# gsRanking(comb.res)

# eaBrowse(comb.res, graph.view=hsa.grn, nr.show=5)

######################################################################
######################################################################
######################################################################
###################################################################### 

# ebrowser( meth=c("ora", "ggea"),
#          exprs=exprs.file, cdat=cdat.file, rdat=rdat.file,
#          org="hsa", gs=hsa.gs, grn=hsa.grn, comb=TRUE, nr.show=5)

######################################################################
###################################################################### 
######################################################################

# configEBrowser(key="OUTDIR.DEFAULT", value="/my/out/dir")

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
