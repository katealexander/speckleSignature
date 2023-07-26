## Makes heatmaps for the top contributing speckle genes for each cancer type

setwd("~/Documents/speckleSignature/Speckle_signature_definition/")

library(pheatmap)

## create new directories to store files if they don't already exist
if (!dir.exists("cBioPortal_signatureHeatmaps_individualCancers")){
  dir.create("cBioPortal_signatureHeatmaps_individualCancers")
}


## load speckle protein rotations
speckles <- read.table("speckleProteinGenes_rotations_allStudiesWithRNA_noNA.txt", sep = "\t", header=T)
row.names(speckles) <- speckles$Gene

for (file in list.files(path = "./speckleProteinGeneExpression")){
  #file <- "pcpg_tcga_pan_can_atlas_2018_speckleExpression.txt" # for testing
  prefix <- strsplit(file, split = "_speckleExpression")[[1]][1]
  expressionDF <- read.table(paste("./speckleProteinGeneExpression/", file, sep=""), sep="\t", header = T, row.names = 1)
  
  ## get top rotation speckle protein genes 
  posGenes <- speckles$Gene.name[speckles[,prefix] > 0.04]
  negGenes <- speckles$Gene.name[speckles[,prefix] < -0.04]
  
  ## subset expression dataframe on the speckle protein gene high contributors
  highContributorGenes <- c(posGenes, negGenes)
  expressionDF <- expressionDF[highContributorGenes,]
  expressionDF <- t(expressionDF) # so that rows are samples, columns are gene
  for (i in 1:ncol(expressionDF)){ # make sure it's all numeric
    expressionDF[,i] <- as.numeric(expressionDF[,i])
  } 
  ## get z scores
  N <- scale(expressionDF, center=TRUE, scale=TRUE)
  N <- as.data.frame(N)
  
  ## column annotations
  posGenes <- as.data.frame(posGenes)
  row.names(posGenes) <- posGenes$posGenes
  colnames(posGenes) <- "Gene"
  negGenes <- as.data.frame(negGenes)
  row.names(negGenes) <- negGenes$negGenes
  colnames(negGenes) <- "Gene"
  posGenes$direction = "pos"
  negGenes$direction = "neg"
  col_annos = rbind(posGenes,negGenes)
  rnames = row.names(col_annos)
  col_annos = as.data.frame(col_annos$direction)
  rownames(col_annos) = rnames
  
  ## color breaks for heatmap
  breaksList = seq(-2, 2, by = .04)
  # pheatmap(N[,1:(length(colnames(N))-1)],treeheight_col = 25, treeheight_row = 25, annotation_row = SpeckleScores, annotation_colors = ann_colors, breaks=breaksList, clustering_method= "ward.D", border_color = NA, cutree_rows = 2, cutree_cols = 2, cex=.8, show_rownames = F, show_colnames = F)
  p = pheatmap(N, 
               treeheight_col = 25, 
               treeheight_row = 25, 
               annotation_col = col_annos, 
               breaks=breaksList, 
               clustering_method= "ward.D", 
               border_color = NA, 
               cutree_rows = 2, 
               cutree_cols = 2, 
               cex=.8, 
               show_rownames = F, 
               show_colnames = F
               )
  heatmapName = paste("cBioPortal_signatureHeatmaps_individualCancers/", prefix, "_heatmap.pdf", sep="")
  pdf(heatmapName, width = 4, height = 4, onefile=FALSE)
  print(p)
  dev.off()
}
