setwd("~/Documents/speckleSignature/Graph_individual_genes/")

library(EnhancedVolcano)

A <- read.table("medianGeneExpression_KIRC_specklepatientGroups.txt", header = T, row.names = 1)
A$log2Ratio <- log2(A$SigI/A$SigII)
published <- read.table("HIF2A_publishedSet.txt")
HIFresponsiveGenes <- read.table("HIF2Atargets_MCF7_786O_combined.txt")
B <- A[A$Gene %in% HIFresponsiveGenes$V1 | A$Gene %in% published$V1,] 
## take only the HIF2A responsive genes that are increasing in both sig I and sig I over normals (but doesn't have to be significant)
B$IvsN <- log2(B$SigI/B$Normal)
B$IIvsN <- log2(B$SigII/B$Normal)
B <- B[(B$IvsN > 0.2 | B$IIvsN > 0.2),] 
B <- B[B$pIvsN < 0.05 | B$pIIvsN < 0.05,] 

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  B$log2Ratio > 0 & B$pIvsII < 0.05, "#CB7245",
  ifelse(B$log2Ratio < -0 & B$pIvsII < 0.05, "#27495C",
         'grey83'))
keyvals[is.na(keyvals)] <- 'grey83'
names(keyvals)[keyvals == '#CB7245'] <- 'Up in Sig I'
names(keyvals)[keyvals == 'grey83'] <- 'n.s.'
names(keyvals)[keyvals == '#27495C'] <- 'Up in Sig II'


p = EnhancedVolcano(B,
                lab = B$GeneSymbol,
                x = 'log2Ratio',
                y = 'pIvsII', 
                boxedLabels = T,
                selectLab = c("ASAP1", "CFLAR", "VEGFA", "DDIT4", "CD70", "PHPT1"),
                #selectLab = c("DGAT1","RB1","GNA13","CFLAR", "SPHK1","NFKBIE","JUN","RPL4","RPL3", "EIF3G","EIF6","E2F1","TFAP4","DDIT4", "ASAP1", "NEDD4", "APOL1", "APOL2", "PHPT1", "EGFR", "KLF7", "PKM", "CD34", "LIMK1", "LDHA", "AURKB", "PLK1", "GALK1", "CLIC4", "BRCA1", "BRCA2", "ABCA1", "ENO1", "ENO2", "VIM", "TGFB1"),
                #title = "Control vs Rai1 KD",
                xlim = c(-2.25, 2.25),
                ylim = c(0, 50),
                pCutoff = 0.00005,
                FCcutoff = 0,
                cutoffLineType = 'blank',
                pointSize = .75,
                labSize = 3,
                colCustom = keyvals,
                colAlpha = 1,
                border = 'full',
                drawConnectors = T,
                arrowheads = FALSE,
                widthConnectors = .5,
                #lengthConnectors = .1,
                gridlines.major = F,
                gridlines.minor = F, 
                borderWidth = 0, 
                labCol = 'grey40',
                colConnectors = 'grey40')

filename = "SigIvsIIbias_volcano.pdf"
pdf(filename, width = 4, height = 5.5, onefile=FALSE)
print(p)
dev.off()
