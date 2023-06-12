## Assigns speckle group to each sample. 1 - Signature I, 2 - tumor adjacent normal, 3 - Signature II
## Has the option to filter out samples based on VHL RNA expression levels. However, I found this didn't make a huge difference, so this option on Line 71 is commented out 

## VHL threshold, see above comment
threshold = 5

## set to reflect name of working directory
setwd("~/Documents/speckleSignature/Graph_individual_genes/")

## set to reflect name of expression file
expressionFile = "~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/KIRC_FPKM_UQ.txt"

library(stringr)
library(dplyr)

# study <- "kirc_tcga_pan_can_atlas_2018"

## get expression levels and put VHL expression into new dataframe
expression <- read.table("~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/KIRC_FPKM_UQ.txt", header = T, row.names = 1, sep = "\t")
VHL <- expression[expression$GENE_SYM == "VHL",2:length(colnames(expression))]
VHL <- t(VHL)
colnames(VHL) <- "VHLexpression"
VHL <- as.data.frame(VHL)
VHL$sampleID <- rownames(VHL)

## get what's needed for speckle score calculation 
Igenes <- read.table("../Calculate_speckle_signature_score/sigI_speckleProteinGenes_22cancers.txt", header = F)
row.names(Igenes) <- Igenes$V1
IIgenes <- read.table("../Calculate_speckle_signature_score/sigII_speckleProteinGenes_22cancers.txt", header = F)
row.names(IIgenes) <- IIgenes$V1
speckles <- read.table("../Calculate_speckle_signature_score/speckleProteinGenes_22cancers.txt", header=T)
row.names(speckles) <- speckles$Gene

##### calculate speckle scores #####
speckleExpression <- expression[expression$GENE_SYM %in% speckles$Gene,]
speckleExpression <- na.omit(speckleExpression)
speckleExpression <- as.data.frame(t(speckleExpression)) # so that rows are samples, columns are gene
colnames(speckleExpression) <- speckleExpression[1,]
speckleExpression <- speckleExpression[2:length(row.names(speckleExpression)),]
for (i in 1:ncol(speckleExpression)){ # make sure it's all numeric
  speckleExpression[,i] <- as.numeric(speckleExpression[,i])
} 
## get z scores
N <- scale(speckleExpression, center=TRUE, scale=TRUE)
N <- as.data.frame(N)
## set speckle signature genes
Igenes$signature = "I"
IIgenes$signature = "II"
## excluding the speckle protein genes not found in dataset and calculating I and II multipliers
Igenes <- Igenes[Igenes$V1 %in% colnames(N),]
IIgenes <- IIgenes[IIgenes$V1 %in% colnames(N),]
Imultiplier = 1.00/length(row.names(Igenes))
IImultiplier = -1.00/length(row.names(IIgenes))
## calculate speckle score, which is the sum of the z scores of Igenes * Imultiplier + sum of the z scores of IIgenes * IImultiplier
N$speckleScore <- rowSums(N[,Igenes$V1]*Imultiplier) + rowSums(N[,IIgenes$V1]*IImultiplier)

## add speckle score to VHL expression dataframe
VHL$speckleScore <- N$speckleScore
VHL$speckleScoreNtile <- ntile(VHL$speckleScore, 4)
VHL$tissueType <- NA
VHL$tissueType[substring(VHL$sampleID, 14,15) == "01"] <- "Primary Tumor"
VHL$tissueType[substring(VHL$sampleID, 14,15) == "11"] <- "Solid Tissue Normal"
VHL$speckleGroup <- NA
VHL$speckleGroup[VHL$speckleScoreNtile == "1"] <- "3"
VHL$speckleGroup[VHL$speckleScoreNtile == "4"] <- "1"
VHL$speckleGroup[VHL$tissueType == "Solid Tissue Normal"] <- "2"
VHL$sampleID <- gsub("\\.", "-", VHL$sampleID)
VHL <- na.omit(VHL)

# ## apply VHL expression thresholds; commented out in this version -- resulsts with or without were similar
# VHL<- rbind(VHL[VHL$tissueType == "Primary Tumor" & VHL$VHLexpression < threshold,], VHL[VHL$tissueType == "Solid Tissue Normal" & VHL$VHLexpression > threshold,])


## write the file
write.table(as.data.frame(VHL[, c(2,6)]), file="KIRC_specklePatientGroups.txt", sep = "\t", quote = F, row.names = F)
 
