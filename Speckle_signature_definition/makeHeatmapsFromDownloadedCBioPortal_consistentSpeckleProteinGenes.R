## Calculates speckle score and adds it to the heatmap

setwd("~/Documents/speckleSignature/Speckle_signature_definition/")

library(pheatmap)
library(viridis)

## create new directories to store files if they don't already exist
if (!dir.exists("cBioPortal_signatureHeatmaps_22cancers")){
  dir.create("cBioPortal_signatureHeatmaps_22cancers")
}

if (!dir.exists("cBioPortal_signatureSpeckleScores_22cancers")){
  dir.create("cBioPortal_signatureSpeckleScores_22cancers")
}

## load speckle protein gene list
speckles <- read.table("speckleProteinGenes_22cancers.txt", header=T)
row.names(speckles) <- speckles$Gene

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
rnames = row.names(col_annos)
col_annos = as.data.frame(col_annos$signature)
rownames(col_annos) = rnames


for (file in list.files(path = "./speckleProteinGeneExpression")){
  prefix <- strsplit(file, split = "_speckleExpression")[[1]][1]
  #file <- "pcpg_tcga_pan_can_atlas_2018_speckleExpression.txt" # for testing
  expressionDF <- read.table(paste("./speckleProteinGeneExpression/", file, sep=""), sep="\t", header = T, row.names = 1)
  expressionDF <- expressionDF[speckles$Gene,]
  expressionDF <- t(expressionDF) # so that rows are samples, columns are gene
  for (i in 1:ncol(expressionDF)){ # make sure it's all numeric
    expressionDF[,i] <- as.numeric(expressionDF[,i])
  } 
  
  ## get z scores
  N <- scale(expressionDF, center=TRUE, scale=TRUE)
  N <- as.data.frame(N)
  
  ## calculate speckle score, which is the sum of the z scores of Igenes * Imultiplier + sum of the z scores of IIgenes * IImultiplier
  N$speckleScore <- rowSums(N[,Igenes$V1]*Imultiplier) + rowSums(N[,IIgenes$V1]*IImultiplier)
  
  ## set row annotations
  SpeckleScores <- N$speckleScore
  SpeckleScores <- as.data.frame(SpeckleScores)
  row.names(SpeckleScores) <- row.names(N)
  
  ## annotation colors
  paletteLength <- 20
  scoreColors <- viridis::viridis(paletteLength)
  ann_colors = list(SpeckleScores = scoreColors)
  
  ## color breaks for heatmap
  breaksList = seq(-2, 2, by = .04)
  # pheatmap(N[,1:(length(colnames(N))-1)],treeheight_col = 25, treeheight_row = 25, annotation_row = SpeckleScores, annotation_colors = ann_colors, breaks=breaksList, clustering_method= "ward.D", border_color = NA, cutree_rows = 2, cutree_cols = 2, cex=.8, show_rownames = F, show_colnames = F)
  p = pheatmap(N[,1:(length(colnames(N))-1)], treeheight_col = 25, treeheight_row = 25, annotation_row = SpeckleScores, annotation_col = col_annos,
               annotation_colors = ann_colors, breaks=breaksList, clustering_method= "ward.D", 
               border_color = NA, cutree_rows = 2, cutree_cols = 2, cex=.8, show_rownames = F, 
               show_colnames = F, annotation_names_row = F)
  heatmapName = paste("cBioPortal_signatureHeatmaps_22cancers/", prefix, "_heatmap.pdf", sep="")
  pdf(heatmapName, width = 4, height = 4, onefile=FALSE)
  print(p)
  dev.off()
  
  write.table(N, sep="\t", file=paste("cBioPortal_signatureSpeckleScores_22cancers/", prefix, "_Zscore_speckleScore.txt", sep=""), quote = F)
  
}
