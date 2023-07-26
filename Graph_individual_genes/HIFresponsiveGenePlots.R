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
published <- read.table("HIF2A_publishedSet.txt")
B <- A[A$Gene %in% HIFresponsiveGenes$V1 | A$Gene %in% published$V1,] 

## take only the HIF2A responsive genes that are increasing in both sig I and sig I over normals
B$IvsN <- log2(B$SigI/B$Normal)
B$IIvsN <- log2(B$SigII/B$Normal)
B <- B[(B$IvsN > 0.3 | B$IIvsN > 0.3),] 
B <- B[B$pIvsN < 0.05 | B$pIIvsN < 0.05,] 

## make heatmap of z-scores
toH <- B[,3:5]
N <- t(scale(t(toH), center=T, scale=T))
p = pheatmap(N, cluster_cols = F, cex = 1, show_rownames = F, clustering_method = "ward.D2", cutree_rows = 5, treeheight_row = 25)
filename = "medianExpression_specklePatientGroups_HIFresponsiveGenes_heatmap.pdf"
pdf(filename, width = 3, height = 4, onefile=FALSE)
print(p)
dev.off()



## cluster compare GO analysis of I-biased and II-biased HIF2A responsive genes, including ns option to change commented
B$log2Ratio <- log2(B$SigI/B$SigII)
Igenes <- B$GeneSymbol[B$log2Ratio > 0.1 & B$pIvsII < 0.05]
IIgenes <- B$GeneSymbol[B$log2Ratio < 0.1 & B$pIvsII < 0.05]
allHIF <- B$GeneSymbol

Igenes.df <- bitr(Igenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
IIgenes.df <- bitr(IIgenes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
allHIF.df <- bitr(allHIF, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db) 

all = enrichGO(gene = allHIF.df$ENTREZID,
               ont  = "BP",
               OrgDb = org.Hs.eg.db,
               pvalueCutoff = 0.01,
               pAdjustMethod = "fdr",
               minGSSize = 5,
               maxGSSize = 300, 
               qvalueCutoff = 0.01,
               readable = FALSE)
all <- clusterProfiler::simplify(all, cutoff = 0.67)
# Convert to data frame and set BP description as rownames
all = as.data.frame(all, row.names = all$Description)

I = enrichGO(gene = Igenes.df$ENTREZID,
                ont  = "BP",
                OrgDb = org.Hs.eg.db,
                pvalueCutoff = 1,
                pAdjustMethod = "fdr",
                minGSSize = 5,
                maxGSSize = 300, 
                qvalueCutoff = 1,
                readable = FALSE)
I = as.data.frame(I, row.names = I$Description)

II = enrichGO(gene = IIgenes.df$ENTREZID,
              ont  = "BP",
              OrgDb = org.Hs.eg.db,
              pvalueCutoff = 1,
              pAdjustMethod = "fdr",
              minGSSize = 5,
              maxGSSize = 300, 
              qvalueCutoff = 1,
              readable = FALSE)
II = as.data.frame(II, row.names = II$Description)

I$sample <- "I"
II$sample <- "II"

ISub <- I[I$Description %in% row.names(all),]
IISub <- II[II$Description %in% row.names(all),]

merged <- merge(ISub, IISub, by = "Description")
sharedTerms <- merged$Description[merged$qvalue.x < 0.05 & merged$qvalue.y < 0.05]
IITerms <- merged$Description[merged$qvalue.x > 0.05 & merged$qvalue.y < 0.05]
ITerms <- merged$Description[merged$qvalue.x < 0.01 & merged$qvalue.y > 0.05]

combo <- rbind(ISub, IISub)
combo$GeneRatio <- combo$Count/as.numeric(str_split(combo$GeneRatio, "/")[[1]][2])
comboSubI <- combo[combo$Description %in% ITerms,]
comboSubShared <- combo[combo$Description %in% sharedTerms,]
comboSubII <- combo[combo$Description %in% IITerms,]


p = comboSubI%>%
  ggplot(aes(x=sample,y=reorder(Description, -log(qvalue)))) +
  geom_point(aes(size= GeneRatio,color= -log(qvalue)) ) +
  #next line add border
  geom_point(aes(size= GeneRatio),shape = 1,colour = "black") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2, limits = c(0, 10), na.value = "red")+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 10),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 10), breaks = seq(0, 0.2, by = .05), limits = c(0,0.2))
filename = "GO_I.pdf"
pdf(filename, width = 6.5, height = 6, onefile=FALSE)
print(p)
dev.off()

p = comboSubII%>%
  ggplot(aes(x=sample,y=reorder(Description, -log(qvalue)))) +
  geom_point(aes(size= GeneRatio,color= -log(qvalue)) ) +
  #next line add border
  geom_point(aes(size= GeneRatio),shape = 1,colour = "black") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2, limits = c(0, 10), na.value = "red")+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 10),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 10), breaks = seq(0, 0.2, by = .05), limits = c(0,0.2))
filename = "GO_II.pdf"
pdf(filename, width = 4.75, height = 2, onefile=FALSE)
print(p)
dev.off()

p = comboSubShared%>%
  ggplot(aes(x=sample,y=reorder(Description, -log(qvalue)))) +
  geom_point(aes(size= GeneRatio,color= -log(qvalue)) ) +
  #next line add border
  geom_point(aes(size= GeneRatio),shape = 1,colour = "black") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2, limits = c(0, 10), na.value = "red")+
  theme_classic() +
  theme(
    legend.position="top",
    legend.title = element_text(face = "bold", size = 8),
    legend.title.align = 2,
    axis.title.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12,colour = 'black', angle = 90),
    axis.text.y = element_text(size = 12,colour = 'black'),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  #if you want to change overall size of dot
  scale_size_continuous(range = c(2, 10), breaks = seq(0, 0.2, by = .05), limits = c(0,0.2))
filename = "GO_shared.pdf"
pdf(filename, width = 4, height = 2, onefile=FALSE)
print(p)
dev.off()

write.table(as.data.frame(Igenes), row.names = F, file = "Igenes.txt", sep="\t", quote = F)
write.table(as.data.frame(IIgenes), row.names = F, file = "IIgenes.txt", sep="\t", quote = F)

