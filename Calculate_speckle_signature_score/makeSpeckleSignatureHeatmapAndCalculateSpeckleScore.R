## Calculates speckle score and adds it to the heatmap

setwd("~/Documents/speckleSignature/Calculate_speckle_signature_score/")

library(pheatmap)
library(viridis)

## load speckle protein gene list
speckles <- read.table("speckleProteinGenes_22cancers.txt", header=T)
row.names(speckles) <- speckles$Gene

## get what's needed for speckle score calculation 
Igenes <- read.table("sigI_speckleProteinGenes_22cancers.txt", header = F)
row.names(Igenes) <- Igenes$V1
IIgenes <- read.table("sigII_speckleProteinGenes_22cancers.txt", header = F)
row.names(IIgenes) <- IIgenes$V1

## get the expression of signature speckle protein genes
expressionDF <- read.table("normalizedCounts_10samples.txt", sep="\t", header = T, row.names = 1)
expressionDF <- expressionDF[speckles$Gene,]
expressionDF <- na.omit(expressionDF)
expressionDF <- as.data.frame(t(expressionDF)) # so that rows are samples, columns are gene

for (i in 1:ncol(expressionDF)){ # make sure it's all numeric
  expressionDF[,i] <- as.numeric(expressionDF[,i])
} 

## get z scores
N <- scale(expressionDF, center=TRUE, scale=TRUE)
N <- as.data.frame(N)

Igenes$signature = "I"
IIgenes$signature = "II"

## excluding the speckle protein genes not found in dataset and calculating I and II multipliers
Igenes <- Igenes[Igenes$V1 %in% colnames(N),]
IIgenes <- IIgenes[IIgenes$V1 %in% colnames(N),]
Imultiplier = 1.00/length(row.names(Igenes))
IImultiplier = -1.00/length(row.names(IIgenes))

## column annotations
col_annos = rbind(Igenes,IIgenes)
rnames = row.names(col_annos)
col_annos = as.data.frame(col_annos$signature)
rownames(col_annos) = rnames
  


## calculate speckle score, which is the sum of the z scores of Igenes * Imultiplier + sum of the z scores of IIgenes * IImultiplier
N$speckleScore <- rowSums(N[,Igenes$V1]*Imultiplier) + rowSums(N[,IIgenes$V1]*IImultiplier)
  
## set row annotations
SpeckleScores <- N$speckleScore
SpeckleScores <- as.data.frame(SpeckleScores)
row.names(SpeckleScores) <- row.names(N)
sorted1 <- order(SpeckleScores$SpeckleScores)
Nordered <- N[sorted1,]

## write table of speckle scores
write.table(as.data.frame(SpeckleScores), sep="\t", file="speckleScores.txt", quote=F)



## annotation colors
paletteLength <- 20
scoreColors <- viridis::viridis(paletteLength)
ann_colors = list(SpeckleScores = scoreColors)

## color breaks for heatmap
breaksList = seq(-2, 2, by = .04)
p = pheatmap(Nordered[,1:(length(colnames(Nordered))-1)], treeheight_col = 25, treeheight_row = 25, annotation_row = SpeckleScores, annotation_col = col_annos,
             annotation_colors = ann_colors, breaks=breaksList, 
             border_color = NA, cutree_rows = 4, cutree_cols = 2, cex=1, show_rownames = T, cluster_rows = F , clustering_method = "complete",
             show_colnames = F, annotation_names_row = F)
file_name = "speckleSignatureHeatmap_PDX.pdf"
pdf(file_name, width = 6.5, height = 3)
print(p)
dev.off()
