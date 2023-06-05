setwd("~/Documents/speckleSignature/Speckle_signature_definition/")

library(pheatmap)

speckle.gene.rotations <- read.table("speckleProteinGenes_rotations_allStudiesWithRNA.txt", sep="\t", header=T, row.names = 1)
speckle.gene.rotations <- as.data.frame(speckle.gene.rotations)
speckle.gene.rotations <- na.omit(speckle.gene.rotations)
write.table(speckle.gene.rotations, sep="\t", file="speckleProteinGenes_rotations_allStudiesWithRNA_noNA.txt", quote = F)
speckle.gene.rotations <- speckle.gene.rotations[15:length(speckle.gene.rotations)]

x <- cor(speckle.gene.rotations, method="pearson")
p = pheatmap(x, cutree_rows = 3, cutree_cols = 3)
pdf("rotationCorrelationsPearson.pdf", width = 8.5, height = 8, onefile=FALSE)
print(p)
dev.off()


