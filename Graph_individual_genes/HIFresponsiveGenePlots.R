setwd("~/Documents/speckleSignature/Graph_individual_genes/")

library(ggplot2)
library(ggpubr)
library(viridis)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

A <- read.table("medianGeneExpression_KIRC_specklepatientGroups.txt", sep="\t", header=T)
HIFresponsiveGenes <- read.table("HIF2Atargets_MCF7_786O_combined.txt")
B <- A[A$Gene %in% HIFresponsiveGenes$V1,] 

## take only the HIF2A responsive genes that are increasing in both sig I and sig I over normals (but doesn't have to be significant)
B$IvsN <- log2(B$SigI/B$Normal)
B$IIvsN <- log2(B$SigII/B$Normal)
B <- B[(B$IvsN > 0 & B$IIvsN > 0),] 

## make heatmap of z-scores
toH <- B[,3:5]
N <- t(scale(t(toH), center=T, scale=T))
p = pheatmap(N, cluster_cols = F, cex = 1, show_rownames = F, clustering_method = "ward.D2", cutree_rows = 5, treeheight_row = 25)
filename = "medianExpression_specklePatientGroups_HIFresponsiveGenes_heatmap.pdf"
pdf(filename, width = 3, height = 4, onefile=FALSE)
print(p)
dev.off()

# ## cluster compare GO analysis of I-biased and II-biased HIF2A responsive genes
# B$log2Ration <- log2(B$SigI/B$SigII)
# Igenes <- B$Gene[B$log2Ratio > 0.15]
# IIgenes <- B$Gene[B$log2Ratio < -0.15]
# unchangingGenes <- B$Gene[B$log2Ratio < 0.15 & B$log2Ratio > -0.15]
# Igenes.df <- bitr(Igenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
# IIgenes.df <- bitr(IIgenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
# geneList <- list(I=Igenes.df$ENTREZID, II=IIgenes.df$ENTREZID)
# y=compareCluster(geneList, fun = "enrichGO", ont = "BP", OrgDb = org.Hs.eg.db)
# Y <- simplify(y, cutoff=0.5,)
# p = dotplot(Y, showCategory = 10)
# filename = "clusterProfiler_IvsIIbias.pdf"
# pdf(filename, width = 6, height = 8, onefile=FALSE)
# print(p)
# dev.off()

## cluster compare GO analysis of I-biased and II-biased HIF2A responsive genes, including ns
B$log2Ratio <- log2(B$SigI/B$SigII)
Igenes <- B$Gene[B$log2Ratio > 0 & B$pIvsII < 0.005]
IIgenes <- B$Gene[B$log2Ratio < 0 & B$pIvsII < 0.005]
unchangingGenes <- B$Gene[B$log2Ratio < 0.15 & B$log2Ratio > -0.15 & B$pIvsII > 0.05]
# unchangingGenes <- B$Gene[B$pIvsII > 0.1]
Igenes.df <- bitr(Igenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
IIgenes.df <- bitr(IIgenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
unchangingGenes.df <- bitr(unchangingGenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
geneList <- list(I=Igenes.df$ENTREZID, ns=unchangingGenes.df$ENTREZID, II=IIgenes.df$ENTREZID)
y=compareCluster(geneList, fun = "enrichGO", ont = "BP", OrgDb = org.Hs.eg.db)
Y <- simplify(y, cutoff=0.5,)
p = dotplot(Y, showCategory =8)
filename = "clusterProfiler_IvsIIbias_withNS.pdf"
pdf(filename, width = 6, height = 7, onefile=FALSE)
print(p)
dev.off()

# cnetplot(Y)


