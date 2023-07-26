## Compares high contributing speckle protein genes between cancer types

setwd("~/Documents/speckleSignature/Speckle_signature_definition/")

library(tidyverse)
library(data.table)
library(pheatmap)
library(VennDiagram)



## load speckle protein rotations
speckles <- read.table("speckleProteinGenes_rotations_allStudiesWithRNA_noNA.txt", sep = "\t", header=T)
row.names(speckles) <- speckles$Gene.name

rotationTable <- speckles[,15:44]
rotationTable <- replace(rotationTable, rotationTable > 0.04, "pos")
rotationTable <- replace(rotationTable, rotationTable < -0.04, "neg")
rotationTable <- replace(rotationTable, rotationTable > -0.04 & rotationTable < 0.04 , "nc")
universe = nrow(rotationTable)

## venns of example
pos1 <- row.names(rotationTable[rotationTable[,"kirc_tcga_pan_can_atlas_2018"] == "pos",])
pos2 <- row.names(rotationTable[rotationTable[,"brca_tcga_pan_can_atlas_2018"] == "pos",])
nOverlap = length(intersect(pos1, pos2))
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(pos1), area2 = length(pos2), cross.area = nOverlap, fill = c("yellow", "#60305d"), cex = 0, alpha = c(0.4,0.9))
filename = "KIRC_BRCA_pos_venny.pdf"
pdf(filename, width = 5, height = 5, onefile=FALSE)
grid.draw(venn.plot)
dev.off()
phyper(nOverlap-1, length(pos1), universe-length(pos1), length(pos2), lower.tail= FALSE)

neg2 <- row.names(rotationTable[rotationTable[,"brca_tcga_pan_can_atlas_2018"] == "neg",])
nOverlap = length(intersect(pos1, neg2))
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(pos1), area2 = length(neg2), cross.area = nOverlap, fill = c("yellow", "Midnight Blue"), cex = 0, alpha = c(0.4,0.9), inverted = T, ext.line.lwd = 0)
filename = "KIRCpos_BRCAneg_venny.pdf"
pdf(filename, width = 5, height = 5, onefile=FALSE)
grid.draw(venn.plot)
dev.off()
phyper(nOverlap, length(pos1), universe-length(pos1), length(neg2))

posOverlaps <- as.data.frame(matrix(ncol = 4, nrow = 0))
negOverlaps <- as.data.frame(matrix(ncol = 4, nrow = 0))


for (study1 in colnames(rotationTable)) {
  pos1 <- row.names(rotationTable[rotationTable[,study1] == "pos",])
  neg1 <- row.names(rotationTable[rotationTable[,study1] == "neg",])
  for (study2 in colnames(rotationTable)) {
    if (study2 != study1){
      pos2 <- row.names(rotationTable[rotationTable[,study2] == "pos",])
      neg2 <- row.names(rotationTable[rotationTable[,study2] == "neg",])
      
      ## positive to positive overlap
      set1 <- paste(study1, "pos", sep = "_")
      set2 <- paste(study2, "pos", sep = "_")
      nOverlap = length(intersect(pos1, pos2))
      p = phyper(nOverlap-1, length(pos1), universe-length(pos1), length(pos2), lower.tail= FALSE, log.p = TRUE)
      posOverlaps = rbind(posOverlaps, c(set1, set2, nOverlap, -p))
      ## negative to negative overlap
      set1 <- paste(study1, "neg", sep = "_")
      set2 <- paste(study2, "neg", sep = "_")
      nOverlap = length(intersect(neg1, neg2))
      p = phyper(nOverlap-1, length(neg1), universe-length(neg1), length(neg2), lower.tail= FALSE, log.p = TRUE)
      negOverlaps = rbind(negOverlaps, c(set1, set2, nOverlap, -p))
      
      ## positive to negative overlap
      
      ## negative to positive ovelap
      
    }
  }
}
colnames(posOverlaps) <- c("set1", "set2", "nOverlap", "phyperInclusion")
colnames(negOverlaps) <- c("set1", "set2", "nOverlap", "phyperInclusion")

## make heatmaps of the p-values
breaksList = seq(0, 10, by = 1)
posTable <- dcast(posOverlaps, set1 ~ set2, value.var = "phyperInclusion")
posTable[upper.tri(posTable)] <- NA
p = pheatmap(sapply(posTable[,2:30], as.numeric), cluster_rows = F, cluster_cols = F, color = inferno(10), na_col = "white", border_color = NA, breaks=breaksList)
filename = "positiveRotation_pvalOfOverlap_heatmap.pdf"
pdf(filename, width = 4, height = 6, onefile=FALSE)
print(p)
dev.off()

negTable <- dcast(negOverlaps, set1 ~ set2, value.var = "phyperInclusion")
negTable[upper.tri(negTable)] <- NA
p = pheatmap(sapply(negTable[,2:30], as.numeric), cluster_rows = F, cluster_cols = F, color = inferno(10), na_col = "white", border_color = NA, breaks = breaksList)
filename = "negativeRotation_pvalOfOverlap_heatmap.pdf"
pdf(filename, width = 4, height = 6, onefile=FALSE)
print(p)
dev.off()
