################################################################################################
################################################################################################

library("enrichR")

################################################################################################
################################################################################################
################################################################################################
################################################################################################

args <- commandArgs(TRUE)

FILE <- args[1] 

name <- basename(FILE)

###### reading the FILE :

file <- read.delim(FILE, sep="\t", header=T, stringsAsFactors=T)

list_genes <- paste("", file$GENE_NAME, "", sep="")

######### FILE <- "analysis.LIMMA.genes.counts.Aph_KH7.only.DEG.and.DOWN.info.all.samples"
      
####################################################################################################
####################################################################################################

library("enrichR")
dbs <- listEnrichrDbs()

### here chosing specific databases :
# dbs <- c("GO_Molecular_Function_2015", 
#         "GO_Cellular_Component_2015", 
#         "GO_Biological_Process_2015")

# enriched <- enrichr(list_genes, dbs)

dbs <- c("GO_Biological_Process_2018",
         "GO_Cellular_Component_2018",
         "GO_Molecular_Function_2018",
         "DSigDB",
         "Genome_Browser_PWMs",
         "TRANSFAC_and_JASPAR_PWMs",
         "ENCODE_TF_ChIP-seq_2014",
         "ENCODE_TF_ChIP-seq_2015",
         "ChEA_2016",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "KEGG_2016",
         "WikiPathways_2016",
         "Reactome_2016",
         "BioCarta_2016",
         "Panther_2016",
         "NCI-Nature_2016", 
         "OMIM_Disease",
         "OMIM_Expanded",
         "MSigDB_Computational",
         "MSigDB_Oncogenic_Signatures",
         "Chromosome_Location")

enriched <- enrichr(list_genes, dbs)

######################################################################################################
######################################################################################################
########## printing the SELECTED databases :  

CATEGORY <- "GO_Biological_Process_2018"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
 
######################################################################################################
######################################################################################################

CATEGORY <- "GO_Cellular_Component_2018"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "GO_Molecular_Function_2018"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "DSigDB"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "Genome_Browser_PWMs"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "TRANSFAC_and_JASPAR_PWMs"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "ENCODE_TF_ChIP-seq_2014"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "ENCODE_TF_ChIP-seq_2015"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "ChEA_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "KEGG_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "WikiPathways_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "Reactome_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "BioCarta_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "Panther_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "NCI-Nature_2016"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "OMIM_Disease"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "OMIM_Expanded"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "MSigDB_Computational"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "MSigDB_Oncogenic_Signatures"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

CATEGORY <- "Chromosome_Location"

results.db <- enriched[[CATEGORY]]

results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

head(results.db.fdr[,1:6])
tail(results.db.fdr[,1:6])
dim(results.db.fdr[,1:6])

write.table(results.db.fdr, 
            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
            sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################

# CATEGORY <- ""

# results.db <- enriched[[CATEGORY]]

# results.db.fdr <- subset(results.db, Adjusted.P.value < 0.05)

# head(results.db.fdr[,1:6])
# tail(results.db.fdr[,1:6])
# dim(results.db.fdr[,1:6])

# write.table(results.db.fdr, 
#            file=paste(basename(FILE), ".enrichR.fdr0.05.", CATEGORY, ".txt", sep=""), 
#            sep="\t", quote = FALSE,
#            row.names = FALSE, col.names = TRUE)

######################################################################################################
#####################################################################################################

#                                 Genome_Browser_PWMs      615        13362
#                            TRANSFAC_and_JASPAR_PWMs      326        27884
#                           Transcription_Factor_PPIs      290         6002
#                                           ChEA_2013      353        47172
#                    Drug_Perturbations_from_GEO_2014      701        47107
#                             ENCODE_TF_ChIP-seq_2014      498        21493
#                                       BioCarta_2013      249         1295
#                                       Reactome_2013       78         3185
#                                   WikiPathways_2013      199         2854
#                Disease_Signatures_from_GEO_up_2014      142        15057
#                                          KEGG_2013      200         4128
#                         TF-LOF_Expression_from_GEO      269        34061
#                                TargetScan_microRNA      222         7504
#                                   PPI_Hub_Proteins      385        16399
#                         GO_Molecular_Function_2015     1136        12753
#                                          GeneSigDB     2139        23726
#                                Chromosome_Location      386        32740
#                                   Human_Gene_Atlas       84        13373
#                                   Mouse_Gene_Atlas       96        19270
#                         GO_Cellular_Component_2015      641        13236
#                         GO_Biological_Process_2015     5192        14264
#                           Human_Phenotype_Ontology     1779         3096
#                    Epigenomics_Roadmap_HM_ChIP-seq      383        22288
#                                           KEA_2013      474         4533
#                  NURSA_Human_Endogenous_Complexome     1796        10231
#                                              CORUM     1658         2741
#                            SILAC_Phosphoproteomics       84         5655
#                    MGI_Mammalian_Phenotype_Level_3       71        10406
#                    MGI_Mammalian_Phenotype_Level_4      476        10493
#                                        Old_CMAP_up     6100        11251
#                                      Old_CMAP_down     6100         8695
#                                       OMIM_Disease       90         1759
#                                      OMIM_Expanded      187         2178
#                                          VirusMINT       85          851
#                               MSigDB_Computational      858        10061
#                        MSigDB_Oncogenic_Signatures      189        11250
#              Disease_Signatures_from_GEO_down_2014      142        15406
#                    Virus_Perturbations_from_GEO_up      323        17711
#                  Virus_Perturbations_from_GEO_down      323        17576
#                      Cancer_Cell_Line_Encyclopedia      967        15797
#                           NCI-60_Cancer_Cell_Lines       93        12232
#        Tissue_Protein_Expression_from_ProteomicsDB      207        13572
#  Tissue_Protein_Expression_from_Human_Proteome_Map       30         6454
#                                   HMDB_Metabolites     3906         3723
#                              Pfam_InterPro_Domains      311         7588
#                         GO_Biological_Process_2013      941         7682
#                         GO_Cellular_Component_2013      205         7324
#                         GO_Molecular_Function_2013      402         8469
#                               Allen_Brain_Atlas_up     2192        13121
#                            ENCODE_TF_ChIP-seq_2015      816        26382
#                  ENCODE_Histone_Modifications_2015      412        29065
#                  Phosphatase_Substrates_from_DEPOD       59          280
#                             Allen_Brain_Atlas_down     2192        13877
#                  ENCODE_Histone_Modifications_2013      109        15852
#                          Achilles_fitness_increase      216         4320
#                          Achilles_fitness_decrease      216         4271
#                       MGI_Mammalian_Phenotype_2013      476        10496
#                                      BioCarta_2015      239         1678
#                                      HumanCyc_2015      125          756
#                                          KEGG_2015      179         3800
#                                    NCI-Nature_2015      209         2541
#                                       Panther_2015      104         1918
#                                  WikiPathways_2015      404         5863
#                                      Reactome_2015     1389         6768
#                                             ESCAPE      315        25651
#                                         HomoloGene       12        19129
#                Disease_Perturbations_from_GEO_down      839        23939
#                  Disease_Perturbations_from_GEO_up      839        23561
#                   Drug_Perturbations_from_GEO_down      906        23877
#                   Genes_Associated_with_NIH_Grants    32876        15886
#                     Drug_Perturbations_from_GEO_up      906        24350
#                                           KEA_2015      428         3102
#              Single_Gene_Perturbations_from_GEO_up     2460        31132
#            Single_Gene_Perturbations_from_GEO_down     2460        30832
#                                          ChEA_2015      395        48230
#                                              dbGaP      345         5613
#                           LINCS_L1000_Chem_Pert_up    33132         9559
#                         LINCS_L1000_Chem_Pert_down    33132         9448
#   GTEx_Tissue_Sample_Gene_Expression_Profiles_down     2918        16725
#     GTEx_Tissue_Sample_Gene_Expression_Profiles_up     2918        19249
#                 Ligand_Perturbations_from_GEO_down      261        15090
#                  Aging_Perturbations_from_GEO_down      286        16129
#                    Aging_Perturbations_from_GEO_up      286        15309
#                   Ligand_Perturbations_from_GEO_up      261        15103
#                   MCF7_Perturbations_from_GEO_down      401        15022
#                     MCF7_Perturbations_from_GEO_up      401        15676
#                Microbe_Perturbations_from_GEO_down      312        15854
#                  Microbe_Perturbations_from_GEO_up      312        15015
#              LINCS_L1000_Ligand_Perturbations_down       96         3788
#                LINCS_L1000_Ligand_Perturbations_up       96         3357
#              LINCS_L1000_Kinase_Perturbations_down     3644        12668
#                LINCS_L1000_Kinase_Perturbations_up     3644        12638
#                                      Reactome_2016     1530         8973
#                                          KEGG_2016      293         7010
#                                  WikiPathways_2016      437         5966
#          ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X      104        15562
#                 Kinase_Perturbations_from_GEO_down      285        17850
#                   Kinase_Perturbations_from_GEO_up      285        17660
#                                      BioCarta_2016      237         1348
#                                     HumanCyc_2016      152          934
#                                   NCI-Nature_2016      209         2541
#                                      Panther_2016      112         2041
#                                        DrugMatrix     7876         5209
#                                         ChEA_2016      645        49238
#                                             huMAP      995         2243
#                                    Jensen_TISSUES     1842        19586
# RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO     1302        22440
#                      MGI_Mammalian_Phenotype_2017     5231         8184
#                               Jensen_COMPARTMENTS     2283        18329
#                                   Jensen_DISEASES     1811        15755
#                                      BioPlex_2017     3915        10271
#                        GO_Cellular_Component_2017      636        10427
#                        GO_Molecular_Function_2017      972        10601
#                        GO_Biological_Process_2017     3166        13822
#                       GO_Cellular_Component_2017b      816         8002
#                       GO_Molecular_Function_2017b     3271        10089
#                       GO_Biological_Process_2017b    10125        13247
#                                    ARCHS4_Tissues      108        21809
#                                 ARCHS4_Cell-lines      125        23601
#                                  ARCHS4_IDG_Coexp      352        20883
#                              ARCHS4_Kinases_Coexp      498        19612
#                                  ARCHS4_TFs_Coexp     1724        25983
#                           SysMyo_Muscle_Gene_Sets     1135        19500
#                                   miRTarBase_2017     3240        14893
#                          TargetScan_microRNA_2017      683        17598
#              Enrichr_Libraries_Most_Popular_Genes      121         5902
#           Enrichr_Submissions_TF-Gene_Coocurrence     1722        12486
#        Data_Acquisition_Method_Most_Popular_Genes       12         1073
#                                            DSigDB     4026        19513
#                        GO_Biological_Process_2018     5103        14433
#                        GO_Cellular_Component_2018      446         8655
#                        GO_Molecular_Function_2018     1151        114

###################################################################################################
####################################################################################################
######### SELECTED DATABASES : ################

# GO_Biological_Process_2018
# GO_Cellular_Component_2018
# GO_Molecular_Function_2018
# DSigDB
# Genome_Browser_PWMs
# TRANSFAC_and_JASPAR_PWMs
# ENCODE_TF_ChIP-seq_2014
# ENCODE_TF_ChIP-seq_2015
# ChEA_2016
# ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
# KEGG_2016
# WikiPathways_2016
# Reactome_2016
# BioCarta_2016   
# Panther_2016
# NCI-Nature_2016 
# OMIM_Disease
# OMIM_Expanded
# MSigDB_Computational
# MSigDB_Oncogenic_Signatures
# Chromosome_Location
