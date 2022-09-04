######################################################################
######################################################################
######################################################################
######################################################################

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(tidyverse)
library(fgsea)

res <- read.delim("analysis.LIMMA.integrating.all.samples.all.genes.in.4.STEPS.on.03oct2019.txt", 
                 header=T, sep="\t", stringsAsFactors=F)

######################################################################
######################################################################

png(paste("figure.barplot", "KEGG", "logFC.COND", "png", sep="."))
barplot(sort(res2$logFC.KH7, decreasing = T))
dev.off()

######################################################################
######################################################################

# deframe() converts two-column data frames to a named vector or list, 
# using the first column as name and the second column as value.

ranks <- deframe(res2)
head(ranks, 20)

###################################################
################################################### GMT_PATHWAYS

pathways.hallmark <- gmtPathways("c2.cp.kegg.v6.2.symbols.gmt")

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100)

######################################################################
######################################################################
######################################################################
######################################################################

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways NES from GSEA") + 
  theme_minimal()

ggsave(paste("figure.barplots", "PATHWAYS.BARPLOTS.KEGG", "logFC.COND", "png", sep="."))

#########################################################################
#########################################################################

# fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100)
# fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100)

head(fgseaRes[order(padj, -abs(NES)), ], n=10)
#                               pathway       pval       padj        ES      NES
# 1:               HALLMARK_P53_PATHWAY 0.01123596 0.03654971 0.7076804 2.792722
# 2:   HALLMARK_TNFA_SIGNALING_VIA_NFKB 0.01136364 0.03654971 0.5758375 2.262639
# 3:                HALLMARK_MYOGENESIS 0.01219512 0.03654971 0.5910913 2.259118

head(fgseaRes[order(padj, -abs(NES)), ], n=10)
sum(fgseaRes[, padj < 0.05])

#######################################################################
#######################################################################
# plotEnrichment(pathway, stats, gseaParam = 1, ticksSize = 0.2)
# Arguments:
# pathway: Gene set to plot.
# stats: Gene-level statistics.
# gseaParam: GSEA parameter.

plotEnrichment(pathways.hallmark[["KEGG_P53_SIGNALING_PATHWAY"]], ranks)
ggsave(paste("figure.GSEA.enrichment", "KEGG", "plot.ENRICHMENT.example.TP53", "logFC.COND", "png", sep="."))

# looking at the genes in these pathways :
# pathways.hallmark %>% 
#   enframe("pathway", "GENE_NAME") %>% 
#   unnest() %>% 
#   inner_join(res, by="GENE_NAME")

#######################################################################
#######################################################################
#######################################################################
####################################################################### GSEA table plot

# The function plotGseaTable allows us to plot a summary figue showing 
# the results for multiple pathways.

topUp <- fgseaRes %>% 
    filter(ES > 0) %>% 
    top_n(10, wt=-padj)

topDown <- fgseaRes %>% 
    filter(ES < 0) %>% 
    top_n(10, wt=-padj)

topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)


png(paste("figure.GSEA.TABLE", "KEGG", "plot.enrichments", "logFC.KH7", "png", sep="."))
plotGseaTable(pathways.hallmark[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)
dev.off()

#######################################################################
#######################################################################
#######################################################################
####################################################################### other sets of PATHWAYS :
#######################################################################
#######################################################################

#fgsea(pathways=gmtPathways("c2.cp.kegg.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
#  as_tibble() %>% 
#  arrange(padj)

#fgsea(pathways=gmtPathways("c3.mir.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
#  as_tibble() %>% 
#  arrange(padj)

#fgsea(pathways=gmtPathways("c5.all.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
#  as_tibble() %>% 
#  arrange(padj)

#######################################################################
#######################################################################
#######################################################################
