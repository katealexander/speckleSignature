setwd("~/Documents/speckleSignature/Calculate_speckle_signature_score/")

library(ggplot2)
library(viridis)
library(pheatmap)
library(dplyr)
library(ggpubr)

## Expression and metadata
expression <- read.table("SRP073253.tsv", sep = "\t", header = T, row.names = 1)
metadata <- read.table("metadata_SRP073253.tsv", sep = "\t", header = T, row.names = 1)

## subcellular localization data downloaded from Human Protein atlas on 02/12/2023
subcellular <- read.table("subcellular_location.tsv", sep = "\t", header =T)
s = "Nuclear speckles"

## get the proteins that are annotated as Enhanced, Supported, or Approved in nuclear speckles
subcellularSpeckle <- subcellular[grep(s, subcellular$Enhanced),]
subcellularSpeckle <- rbind(subcellularSpeckle, subcellular[grep(s, subcellular$Supported),])
subcellularSpeckle <- rbind(subcellularSpeckle, subcellular[grep(s, subcellular$Approved),])
row.names(subcellularSpeckle) <- subcellularSpeckle$Gene.name

## get what's needed for speckle score calculation 
Igenes <- read.table("sigI_speckleProteinGenes_22cancers.txt", header = F)
row.names(Igenes) <- Igenes$V1
IIgenes <- read.table("sigII_speckleProteinGenes_22cancers.txt", header = F)
row.names(IIgenes) <- IIgenes$V1
Imultiplier = 1.00/length(row.names(Igenes))
IImultiplier = -1.00/length(row.names(IIgenes))

## column annotations
Igenes$signature = "I"
IIgenes$signature = "II"
col_annos = rbind(Igenes,IIgenes)
geneSymbol = row.names(col_annos)
col_annos = as.data.frame(col_annos$signature)
col_annos$geneSymbol = geneSymbol
row.names(col_annos) = subcellularSpeckle$Gene[match(col_annos$geneSymbol, subcellularSpeckle$Gene.name)]

## get speckle protein gene expression
expressionSpeckleSignature <- expression[row.names(col_annos),]
expressionSpeckleSignature <- t(expressionSpeckleSignature) # so that rows are samples, columns are gene
expressionSpeckleSignature <- as.data.frame(expressionSpeckleSignature)

## get z scores
N <- scale(expressionSpeckleSignature, center=TRUE, scale=TRUE)
N <- as.data.frame(N)

## calculate speckle score, which is the sum of the z scores of Igenes * Imultiplier + sum of the z scores of IIgenes * IImultiplier
N$speckleScore <- rowSums(N[,row.names(col_annos)[col_annos$`col_annos$signature`=="I"]]*Imultiplier) + rowSums(N[,row.names(col_annos)[col_annos$`col_annos$signature`=="II"]]*IImultiplier)
N$drug_response <- metadata$drug_response[match(row.names(N), row.names(metadata))]
N$tissue <- metadata$refinebio_specimen_part[match(row.names(N), row.names(metadata))]

## set row annotations
row_annos <- N$speckleScore
row_annos <- as.data.frame(row_annos)
colnames(row_annos) <- "speckleScore"
row.names(row_annos) <- row.names(N)
row_annos$drug_response <- N$drug_response
row_annos$tissue <- N$tissue

metadata$speckleScore <- N$speckleScore[match(row.names(N), row.names(metadata))]
write.table(metadata, sep = "\t", file = "metadata_SRP073253_withSpeckleScore.txt", quote = F)




