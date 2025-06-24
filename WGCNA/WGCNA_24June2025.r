# Datasets : 

# https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html
# https://fuzzyatelin.github.io/bioanth-stats/module-25/module-25.html

# https://edu.sib.swiss/pluginfile.php/158/course/section/65/_01_SIB2016_wgcna.pdf
# https://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html
# Horvath S (2011) Weighted Network Analysis. Applications in Genomics and Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8

library(tidyverse)
library(magrittr)      
library(WGCNA)        
library(flashClust)

# flashClust is an R package that provides a faster implementation of hierarchical clustering (hclust) for large datasets, 
# particularly useful in applications like WGCNA.

setwd("/home/tanasa/WGCNA")
list.files()

# the article : https://pmc.ncbi.nlm.nih.gov/articles/PMC2631488/

# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61333/suppl/GSE61333_ligule_count.txt.gz
# gunzip GSE61333_ligule_count.txt.gz
# data <- readr::read_delim("data/GSE61333_ligule_count.txt", delim = "\t")

# https://rpubs.com/natmurad/WGCNA

# working with :
# wget 'https://raw.githubusercontent.com/fuzzyatelin/fuzzyatelin.github.io/master/bioanth-stats/module-F21-Group1/FemaleLiver-Data/LiverFemale3600.csv'
# https://raw.githubusercontent.com/fuzzyatelin/fuzzyatelin.github.io/master/bioanth-stats/module-F21-Group1/FemaleLiver-Data/LiverMale3600.csv'
# https://raw.githubusercontent.com/fuzzyatelin/fuzzyatelin.github.io/master/bioanth-stats/module-F21-Group1/FemaleLiver-Data/ClinicalTraits.csv

# GPT4 : an overview of WCGNA

# ⚙️ WGCNA Algorithm: Step-by-Step

# Input :
#    A gene expresion matrix: genes × samples
#        Rows: (or transcripts)
#        Columns: samples (from conditions, timepoints, individuals, etc.)
#        Matrix should be normalized (e.g., log2(CPM + 1), vst, rlog)

# 🔢 Step 1: Compute pairwise correlations

# Input: matrix G (genes x samples)
# Output: symmetric correlation matrix S (genes x genes)
#
#    Compute Pearson correlation between all gene pairs:
#    sij=cor(xi,xj)
#    sij​=cor(xi​,xj​)

#    Matrix S=[sij]S=[sij​] is the similarity matrix

# 🔁 Step 2: Construct adjacency matrix using soft thresholding

# Input: similarity matrix S
# Output: adjacency matrix A (genes x genes)

#    Apply a power adjacency function to emphasize strong correlations:
#    aij=∣sij∣β
#    aij​=∣sij​∣β

#  where β i a soft-thresholding power, chosen so the resulting network approximates a scale-free topology 
# (i.e., few hubs with many connections).

#    this step produces the adjacency matrix AA, encoding the connection strength between each gene pair.

# 🔄 Step 3: Compute Topological Overlap Matrix (TOM)

# Input: adjacency matrix A
# Output: TOM matrix T (genes x genes)

#    TOM quantifies shared neighbors between gene pairs and provides a more biologically robust similarity measure:
#    TOMij=lij+aijmin⁡(ki,kj)+1−aij
#    TOMij​=min(ki​,kj​)+1−aij​lij​+aij​​

#    where:
#        lij=∑uaiu⋅ajulij​=∑u​aiu​⋅aju​: number of shared neighbors
#        ki=∑uaiuki​=∑u​aiu​: connectivity of node ii

#    Result: a smoother, more biologically meaningful similarity matrix

# 🌳 Step 4: Hierarchical clustering of genes

# Input: dissimilarity = 1 - TOM
# Output: dendrogram of genes

#  Use 1−TOM1 : TOM as a distance matrix

#  Apply hierarchical clustering (average linkage)
#  The result is a dendrogram representing gene co-expression structure

# 🎨 Step 5: Identify modules via dynamic tree cut

# Input: gene dendrogram
# Output: module labels (e.g., colors or cluster IDs)

#    Cut the dendrogram using dynamic tree cut:
#        Identifies branches representing distinct gene co-expression modules
#        Parameters control minimum module size, deep split level, etc.
#    Each gene is assigned a module label (color or number)

# 📈 Step 6: Compute module eigengenes

# Input: module labels, original expression matrix
# Output: module eigengenes (1 per module)

#    For each module:

#        Calculate the module eigengene (ME) — first principal component of the expression profiles of genes in the module:
#        MEk=PC1(Gk)
#        MEk​=PC1​(Gk​)

#    The ME summarizes the module’s overall expression pattern across samples

# 🧪 Step 7: Relate modules to traits

# Input: module eigengenes, sample traits
# Output: module–trait correlation matrix

#    Correlate each module eigengene with external sample traits (e.g., clinical, environmental, experimental)
#    cor(MEk,Trait)
#    cor(MEk​,Trait)

#    Produces:
#        Correlation matrix
#        P-value matrix
#        Heatmap visualization

# 🧩 Step 8: Identify key genes (hub genes)

#    Module Membership (MM): correlation of a gene with its module eigengene
#    MMik = cor(xi,MEk)
#    MMik​ = cor(xi​,MEk​)

#    Gene Significance (GS): correlation between a gene’s expression and the trait
#    GSi=cor(xi,Trait)
#    GSi​=cor(xi​,Trait)

#    Genes with high MM + high GS are considered hub genes in biologically relevant modules

# 🧮 Summary of Core Matrices :

# Matrix	Shape	Meaning
# Expression matrix	         genes × samples	Normalized input
# Similarity (S)	         genes × genes	Pairwise correlation
# Adjacency (A)	             genes × genes	Soft-thresholded similarity
# TOM	                     genes × genes	Shared neighborhood similarity
# Dissimilarity	             genes × genes	1−TOMij1−TOMij​
# Module eigengenes	samples × modules	PC1 of each module

# The following script WGCNA.r performs a Weighted Gene Co-expression Network Analysis (WGCNA) on gene expression data 
# from female mouse livers, aiming to identify gene modules and relate them to clinical traits.

# Google Gemini  : Here's a breakdown of the script's key steps and what each part accomplishes:

# 1. Data Loading and Preprocessing

#    Loads gene expression data from LiverFemale3600.csv.
#    Removes non-expression columns and transposes the data so that rows are samples and columns are genes.
#    Renames columns with gene IDs.

# 2. Outlier Detection

#    Uses goodSamplesGenes() to identify and remove genes and samples with too many missing values or zero variance.
#    Performs hierarchical clustering on samples to visually identify and remove outlier samples (e.g., "F2_221").

# 3. Network Construction: Choosing the Soft-Thresholding Power (β)

#    Explains that WGCNA uses weighted networks where connection strength is determined by a power function of the Pearson correlation coefficient.
#    The pickSoftThreshold() function is used to calculate a suitable soft-thresholding power (β). This parameter is crucial for constructing a scale-free network, where a few genes (hubs) have many connections and most genes have few.
#    The script plots the scale-free topology model fit (R2) and mean connectivity against different soft threshold powers to help in selecting the optimal β (aiming for R2 > 0.8 and minimizing connectivity loss). In this case, β=6 is chosen.
#    The adjacency() function then calculates the adjacency matrix using the chosen soft power.

# 4. Module Construction

#    Topological Overlap Matrix (TOM): The adjacency matrix is transformed into a Topological Overlap Matrix (TOM) and then into a TOM-based dissimilarity matrix (1 - TOM). TOM is preferred because it considers shared neighbors, leading to more biologically distinct modules.
#    Hierarchical Clustering: Hierarchical clustering is applied to the TOM dissimilarity matrix to create a gene dendrogram.
#    Dynamic Tree Cut: The cutreeDynamic() function is used to cut the dendrogram and identify gene modules, assigning each gene a module label (represented by a color). A minimum module size of 30 is recommended.
#    Module Eigengenes (MEs): moduleEigengenes() calculates the module eigengene for each module. An ME is the first principal component of the expression profiles of genes within a module, representing the module's overall expression pattern.

# 5. Module Merging

#    Modules with similar expression profiles (based on the correlation of their eigengenes) are merged using mergeCloseModules(). A height cut-off of 0.25 (corresponding to a correlation of 0.75) is used for merging.

# 6. Relating Modules to External Traits

#    Loads clinical trait data from ClinicalTraits.csv.
#    Matches and aligns the sample names between the expression data and the clinical trait data.
#    Calculates the Pearson correlation and corresponding p-values between each merged module eigengene and the available clinical traits.
#    Visualizes these relationships using a labeledHeatmap(), showing strong positive correlations in red and strong negative correlations in blue.

# 7. Target Gene Identification (Hub Genes)

#    Gene Significance (GS): Calculates the correlation between each gene's expression and a specific trait (e.g., "weight_g").
#    Module Membership (MM): Calculates the correlation of each gene with its assigned module eigengene.
#    Genes with high module membership and high gene significance for a trait are considered potential hub genes within biologically relevant modules.
#    Scatter plots are generated to visualize the relationship between module membership and gene significance for specific modules (e.g., "magenta" and "purple").
#    plotEigengeneNetworks() is used to visualize the relationships among module eigengenes and traits.

d = "LiverFemale3600.csv"

# This .csv file contains quantified and normalized gene expression data from livers of several F2 female mice.

# GPT4: 

# The LiverFemale dataset in WGCNA is a built-in example dataset that contains female mouse liver gene expression data, 
# used primarily for teaching and demonstrating the workflow of WGCNA (Weighted Gene Co-expression Network Analysis).

# 🧬 What it contains:
# It is a microarray gene expression dataset derived from mouse liver tissue.

# Specifically, it contains expression levels from female mice, often across a genetically diverse panel 
# (e.g., from the BXD recombinant inbred strains).

# Often paired with a traitData dataset containing phenotypic traits (like body weight, fat percentage, etc.) for correlation analyses.

# ✅ What it is good for:
# Practicing the WGCNA pipeline, including:
# Data preprocessing
# Network construction
# Module detection
# Relating modules to external traits
# Hub gene analysis
# Understanding module-trait relationships, since the dataset includes phenotype traits like:
# Body weight
# Fat percentage
# Cholesterol levels
# and others

# 📚 How to load it:
# If using the official WGCNA tutorials:
# library(WGCNA)
# data(LiverFemale)

liver.data <- read.csv(file = d, stringsAsFactors = FALSE, header = TRUE)
head(liver.data, 2)
tail(liver.data, 2)

expression.data <- liver.data[,-c(1:8)] # removing variables not holding expression data
head(expression.data, 2)
tail(expression.data, 2)

expression.data <- as.data.frame(t(expression.data)) # transforming the data.frame so columns now represent genes and rows represent samples
names(expression.data) <- liver.data$substanceBXH    # renaming the columns so we don't lose the gene IDs

# Identifying Outlier Genes

# The WGCNA package has a built in function to identify outlier genes called goodSampleGenes(). 
# The function checks the data and returns a list object of samples and genes that pass its filtering criteria

head(expression.data, 2)
tail(expression.data, 2)

gsg <- goodSamplesGenes(expression.data)
summary(gsg)
gsg$allOK

# Should we have to remove samples or genes :

if (!gsg$allOK)
{

# Identifies and prints outlier genes
if (sum(!gsg$goodGenes)>0) 
printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); 

# Identifies and prints oulier samples
# Removes the offending genes and samples from the data

if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); 
    
expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] 
}

# Identifying outlier samples

# You can identify outlier samples by using hierarchical clustering. 

sampleTree <- hclust(dist(expression.data), method = "average") # Clustering samples based on distance 

# Setting the graphical parameters
par(cex = 0.6)          # cex controls the scaling of text elements (like axis labels, titles, etc.).
par(mar = c(0,4,2,0))   # c(bottom, left, top, right) — so this line means:

# Plotting the cluster den drogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

# Here you can see sample F2_221 seems to be distant from all of the other samples indicating it is likely an outlier sample.

sampleTree

# Here you can see sample F2_221 seems to be distant from all of the other samples indicating it is likely an outlier sample.
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10) # returns numerical vectors

# Remove outliers
expression.data <- expression.data[cut.sampleTree==1, ]
dim(expression.data)

print("Network Construction")

print("Pairwise Gene Co-expression similarity")

# The WGCNA the similarity measurement for each pair of genes (gene i and gene j) is denoted by their Pearson correlation coefficient.

# Unsigned Networks
# For unsigned networks you take the absolute value of the correlation. 
# If the network is unsigned, you do not know 
# if the gene is up or down regulated in the sample trait, just that its expression is significantly different.

# Signed Networks

print("Adjacency: Pairwise connection")

# The next step after calculating the similarity measurement for each gene pair is to translate that similarity measurement 
# into the gene pairs adjacency to generate an adjacency matrix. 
# Adjacency is the assignment of a connection strength based on the co-expression similarity measurement (Pearson correlation coefficient). 
# Nodes are considered connected if they have a significant pairwise correlation association.

# In unweighted networks, the adjacency matrix indicates whether or not a pair of nodes are connected in a binary fashion.

# Weighted Networks
# In weighted networks the adjacency/connection is not binary and therefore can also distinguish the strength of connection.
# Weighted networks utilize a power function based on a soft threshold parameter 𝛽

# To determine the 𝛽 parameter in a weighted network analysis we try to maximize a model fit (𝑅2 value under a linear regression analysis) 
# under a scale free topology model, while minimizing the number of connections lost when fitting the model 
# (maintaining a high mean number of connections). 
# As 𝑅2 values approach 1, we usually see networks with very few connections. Usually this happy medium occurs at 𝑅2 of > .8.

# The scale-free topology is used because it is based on the idea that the probability that a node (gene) is connected with k other nodes 
# (genes) decays as a power law: 𝑝(𝑘)∼𝑘−𝛾

# In network theory, a scale-free topology refers to a type of network where:
# 🔑 Some nodes (called hubs) have many more connections than others, and
# 📈 The node degree distribution follows a power law : P(k) ∝ k⁻ᵞ, where P(k) is the probability a node has k connections.

# The pickSoftThreshold() function calculates multiple networks all based on different 𝛽 values and returns a data frame 
# with the 𝑅2 values for the networks scale-free topology model fit as well as the mean connectivity measures.

spt <- pickSoftThreshold(expression.data) 
str(spt)

# You can then plot this data frame to better visualize what 𝛽 value you should choose.
# REMINDER : we should be maximizing the 𝑅2 value and minimizing mean connectivity.

# Plot the 𝑅2 values as a function of the soft thresholds.

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2", 
     type="n",
     main = paste("Scale independence"))

text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

# Plot mean connectivity as a function of soft thresholds

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], 
     spt$fitIndices[,5],
     xlab="Soft Threshold (power)", 
     ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

# You can determine the soft power threshold should be set to 6 as it is the spt that retains the highest mean connectivity 
# while reaching an 𝑅2 value above 0.80.



print("Calling the Adjacency Function")

# Now that you have the soft threshold power determined you can call on the adjacency() function of the WGCNA package.
# This function calculates the similarity measurement and transforms the similarity by the adjacency function 
# and generates a weighted network adjacency matrix.

softPower <- 6
adjacencyResult <- adjacency(expression.data, power = softPower)

head(adjacencyResult, 2)
tail(adjacencyResult, 2)



print("Module Construction : Defining Dissimilarity")

# Once the network is constructed, you can begin to extract some meaningful relationships. 
# You can use hierarchical clustering yet again to cluster the network into modules.

# NOTE: A module is a group of gene profiles that are highly correlated, or have a high topological overlap.

# In order to utilize the clustering functions in R you must transform the adjacency matrix into measures of gene dissimilarity 
# (distance of a gene from every other gene in the system).

# NOTE: This is due to the fact that dissimilarity is used in traditional cluster analyses.

print("Topological Overlap Matrix")

# The TOM-based dissimilarity is preferentially used over dissimilarity based on correlation coefficients because 
# TOM-based dissimilarity generates more distinct modules.

# TOMsimilarity # Calculation of the topological overlap matrix, and the corresponding dissimilarity, from a given adjacency matrix. 

TOM <- TOMsimilarity(adjacencyResult)
head(TOM, 2)
tail(TOM, 2)

# To convert this matrix into a dissimilarity matrix you can subtract the TOM object from 1.

TOM.dissimilarity <- 1 - TOM

print("Hierarchical Clustering Analysis")

# The dissimilarity/distance measures are then clustered using linkage hierarchical clustering and a dendrogram (cluster tree) of genes 
# is constructed.

# creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
geneTree

# plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", 
     main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, 
     hang = 0.04)

# To identify modules from this gene dendrogram, you can use the cutreeDynamic() function : it will allow you to set a minimum cluster size. 
# For genomics data it is more beneficial to set minimum module sizes relatively high as you are working with high loads of data. 
# The authors of WGCNA recommend to start at a minClusterSize = 30.

Modules <- cutreeDynamic(dendro = geneTree, 
                         distM = TOM.dissimilarity, 
                         deepSplit = 2, 
                         pamRespectsDendro = FALSE, 
                         minClusterSize = 30)

head(Modules)
tail(Modules)

table(Modules) # it returns a table of the counts of factor levels in an object.
               # i.e.  how many genes are assigned to each created module. 

ModuleColors <- labels2colors(Modules) # it assigns each module number a color
table(ModuleColors)                    # it returns the counts for each color (aka the number of genes within each module)

# it plots the gene dendrogram with the module colors

plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram")

print("Module Eigengene Identification")

# A ME (Module Eigengene) is the standardized gene expression profile for a given module.
# To identify the Module Eigengene you can call on the expression data into the moduleEigengenes() function.

MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs, 2)
tail(MEs, 2)

print("Module Merging")

# To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations 
# and merge the modules that have similar expression profiles.
# REMINDER: An eigengene is the gene whose expression is representative of the the majority of genes expressed within a module.

ME.dissimilarity = 1 - cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

METree = hclust(as.dist(ME.dissimilarity), method = "average") # Clustering eigengenes 
par(mar = c(0,4,2,0))                                          # seting margin sizes
par(cex = 0.6)                                                 # scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

# This figure shows all of the modules which are more than 75% similar. 
# For example you can see that MEcyan and MEpurple are more than 75% similar. 
# Now you can merge these two modules, and others like them using the mergeCloseModules().

# Merge modules whose dissimilarity is below the cutoff !
# corFnc=“pearson”; power=6; min. module size=30

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
mergedColors = merge$colors

# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# a dendrogram which shows both the original AND merged module colors

plotDendroAndColors(geneTree, 
                    cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")



print("Matching the External Traits")

# Compute correlations: each module eigengene to each trait variable cor(MEs, traitDat)

traitData <- read.csv("ClinicalTraits.csv", header = TRUE, stringsAsFactors = FALSE)
head(traitData, 2)

allTraits <- traitData[, -c(31, 16)]   # removing notes and comments sections 
allTraits <- allTraits[, c(2, 11:36) ] # pulling out only continuous traits 
head(allTraits, 2)
tail(allTraits, 2)

# you must match the trait data to the expression data by the sample number

# 1. Get sample names from expression matrix
Samples <- rownames(expression.data)

# 2. Clean sample names to avoid case/space mismatch
Samples_clean <- trimws(tolower(Samples))
Traits_clean <- trimws(tolower(allTraits$Mice))

# 3. Match cleaned sample names
traitRows <- match(Samples_clean, Traits_clean)

# 4. Identify unmatched samples (if any)
if (any(is.na(traitRows))) {
  cat("⚠️ The following samples were not matched in allTraits$Mice:\n")
  print(Samples[is.na(traitRows)])

  # Optionally drop unmatched samples from expression matrix
  validSamples <- Samples[!is.na(traitRows)]
  expression.data <- expression.data[validSamples, ]
  Samples <- rownames(expression.data)

  # Re-run match on filtered samples
  Samples_clean <- trimws(tolower(Samples))
  traitRows <- match(Samples_clean, Traits_clean)
}

# 5. Extract traits and set row names to match expression data
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

# 6. Final alignment sanity check
stopifnot(identical(rownames(datTraits), rownames(expression.data)))

# 7. Optional: view result
head(datTraits, 2)
tail(datTraits, 2)
dim(datTraits)
dim(expression.data)

identical(rownames(datTraits), rownames(expression.data))

identical(rownames(mergedMEs), rownames(datTraits))

# Fix row order of datTraits
datTraits <- datTraits[match(rownames(mergedMEs), rownames(datTraits)), ]
stopifnot(identical(rownames(mergedMEs), rownames(datTraits)))



# Table of module – trait correlations
# • Identify modules highly correlated to traits of interest
# • Identify traits highly correlated to multiple modules

print("Module-Trait associations")

# Define numbers of genes and samples
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)

nGenes
nSamples

module.trait.correlation = cor(mergedMEs, datTraits, use = "p")            # p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) # calculate the p-value associated with the correlation

head(module.trait.correlation, 2)
head(module.trait.Pvalue, 2)

# Will display correlations and their p-values

textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))

# Display the correlation values within a heatmap plot

labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Each row corresponds to a module eigengene, and the columns correspond to a trait. Each cell contains a p-value and correlation. 
# Those with strong positive correlations are shaded a darker red while those with stronger negative correlations become more blue.



print("Target Gene Identification")

# You can use the gene significance along with the genes intramodular connectivity to identify potential target genes associated 
# with a particular trait of interest. For this analysis weight will be the clinical trait.

# Connectivity - how connected a speficic node is in the network (how many nodes have high correlation with that node). 
# High connectivity indicates a hub gene (central to many nodes). 

# Whole Network connectivity - a measure for how well the node is connected throughout the entire system 
# Intramodular connectivity - a measure for how well the node is connected within its assigned module. 
# Also an indicator for how well that node belongs to its module. This is also known as module membership.

# Potential driver genes
# Strategy:
# Identify those genes within a module that are
# 1) Highly connected within the module (hub genes)
# AND
# 2) Most strongly correlated with a clinical/phenotypical trait of interest

# How to detect hub genes inside a module?

# The straightforward way: gene(s) with highest intramodular connectivity (= sum of in-module edge weights)
# Alternative way proposed in WGCNA: gene(s) with highest module membership

# Module Membership of Genes

# Module membership: Correlation of a gene to a module eigengene
# • Genes with high module membership are good representatives of the overall expression profile in the module
# • Genes with high module membership tend to be “hub” genes in the module (high intramodule connectivity)
# • A gene can have high membership in several modules (not just the one to which it is assigned)

# Driver genes : most correlated with both module and trait (picked by visual inspection)

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

# Define variable weight containing the weight column of datTrait

weight = as.data.frame(datTraits$weight_g)
names(weight) = "weight"

modNames = substring(names(mergedMEs), 3) #e xtract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
head(MMPvalue, 2)

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")
head(GSPvalue, 2)

# Using the gene significance you can identify genes that have a high significance for weight. 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

# Goal:
# To visualize the relationship between:
# Module Membership (MM): How strongly each gene belongs to a module (e.g., magenta)
# Gene Significance (GS): How correlated each gene is with a trait (e.g., body weight)
# This helps identify hub genes: genes that are both central to a module and strongly related to a phenotype.

# As an example, you can look at the magenta module as it has the highest significant association with weight (.59).
# Plot a scatter plot of gene significance vs. module membership in the magenta module.

par(mar=c(1,1,1,1))
module = "magenta"
column = match(module, modNames)
moduleGenes = mergedColors==module # Creates a logical vector identifying which genes belong to the magenta module.

verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# The magenta gene significance and module membership have a positive correlation of .49 with a very significant p-value. 
# This indicates that the genes that are highly significantly associated with the trait (high gene significance) are also the genes 
# that are the most connected within their module (high module membership). 
# Therefore genes in the magenta module could be potential target genes when looking at body weight.

library(ggplot2)

# another way to display the data 
df_plot <- data.frame(
  ModuleMembership = abs(geneModuleMembership[moduleGenes, column]),
  GeneSignificance = abs(geneTraitSignificance[moduleGenes, 1])
)

ggplot(df_plot, aes(x = ModuleMembership, y = GeneSignificance)) +
  geom_point(color = module, alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  labs(
    title = paste("Module membership vs. gene significance\n"),
    subtitle = paste("Module:", module),
    x = paste("Module Membership in", module, "module"),
    y = "Gene significance for body weight"
  ) +
  theme_minimal(base_size = 14)


# Considering another module - "purple"
# Does the module you selected contain potential candidate genes ?

par(mar=c(1,1,1,1))
module = "purple"
column = match(module, modNames)
moduleGenes = mergedColors==module

verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                  abs(geneTraitSignificance[moduleGenes,1]),
                  xlab = paste("Module Membership in", module, "module"),
                  ylab = "Gene significance for body weight",
                  main = paste("Module membership vs. gene significance\n"),
                  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# ggplot2 version of the figure :

library(ggplot2)

# 1. Define module and match relevant column
module <- "purple"
column <- match(module, modNames)
moduleGenes <- mergedColors == module

# 2. Create a data frame with absolute values
df_plot <- data.frame(
  ModuleMembership = abs(geneModuleMembership[moduleGenes, column]),
  GeneSignificance = abs(geneTraitSignificance[moduleGenes, 1])
)

# 3. Compute correlation and p-value
cor_val <- cor(df_plot$ModuleMembership, df_plot$GeneSignificance)
p_val <- cor.test(df_plot$ModuleMembership, df_plot$GeneSignificance)$p.value

# 4. Create ggplot
ggplot(df_plot, aes(x = ModuleMembership, y = GeneSignificance)) +
  geom_point(color = module, size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(
    title = "Module membership vs. gene significance",
    subtitle = paste("Module:", module),
    x = paste("Module Membership in", module, "module"),
    y = "Gene significance for body weight"
  ) +
  annotate(
    "text",
    x = 0.1, y = max(df_plot$GeneSignificance),
    label = sprintf("cor = %.2f, p = %.1e", cor_val, p_val),
    hjust = 0, size = 5
  ) +
  theme_minimal(base_size = 14)




print("Network Visualization of Eigengenes")
# It is also possible to study the relationship among found modules. 
# One way to do this is to quantify module similarity (adjacency) by calculating the pairwise correlation of representative eigengenes.

# You are able to generate a summary plot of the eigengene correlations using the function plotEigengeneNetworks(). 
# To visualize how the trait fits into the eigengene network you can bind the trait data to the module eigengenes 
# and argue the data frame to plotEigengeneNetworks() which will display the relationship on a heatmap.

# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))

# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(5,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle = 90)

# With this heatmap you can identify groups of correlated eigengenes called meta modules. 
# Modules with mutual correlations stronger than their correlation with the specified clinical trait would be grouped into a meta module.

# Another way to show the data :

# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", 
                      marHeatmap = c(5,5,2,2),
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)



# Create a dendrogram and heatmap using two clinical traits. For example, weight AND Glucose.

# Isolate weight and Glucose from the clinical traits

weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
Glucose = as.data.frame(datTraits$Glucose)
names(Glucose) = "Glucose"

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight,Glucose))

# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(5,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle = 90)



# Quality checks on modules :

# Connectivity
# mean intra-module connectivity
# mean ratio of intra-module / total connectivity
# Trait correlations
# strong correlation between module eigengenes and traits of interest
# strong correlation between gene module membership and gene-trait corr.
# Functional enrichment
# many functionally related genes in the same module : GO / pathway analysis

# Assessing modules by intramodular connectivity
# • Simple ranking of modules (highest to lowest mean connectivity)
# • Statistical analysis : obtain p-values, e.g. via bootstrap




