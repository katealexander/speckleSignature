setwd("~/Documents/speckleSignature/Calculate_speckle_signature_score/")

library(ggplot2)
library(viridis)
library(pheatmap)
library(dplyr)
library(ggpubr)
library(survival)
library(survminer)


## Expression and metadata
expression <- read.table("~/Desktop/HIF and Speckles/Braun2020_Checkmate_RNAseq_SpeckleSignature/NIHMS1611472-supplement-Table_S4__RNA_expression__normalized_expression_matrix_.txt", sep = "\t", header = T) ## this file was too big for github, so will need to be downloaded freshly form Braun et al., 2020
row.names(expression) <- make.unique(expression$gene_name)
metadata <- read.table("NIHMS1611472-supplement-Table_S1__Clinical_and_immune_phenotype_data_for_the_CheckMate_cohorts.txt", sep = "\t", header = T, row.names = 1)

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
Igenes <- Igenes[intersect(row.names(expression), Igenes$V1),]
Igenes <- as.data.frame(Igenes)
row.names(Igenes) <- Igenes$Igenes
IIgenes <- read.table("sigII_speckleProteinGenes_22cancers.txt", header = F)
row.names(IIgenes) <- IIgenes$V1
IIgenes <- IIgenes[intersect(row.names(expression), row.names(IIgenes)),]
IIgenes <- as.data.frame(IIgenes)
Imultiplier = 1.00/length(row.names(Igenes))
IImultiplier = -1.00/length(row.names(IIgenes))

## column annotations
Igenes$signature = "I"
names(Igenes) = c("geneSymbol", "signature")
IIgenes$signature = "II"
names(IIgenes) = c("geneSymbol", "signature")
col_annos = rbind(Igenes,IIgenes)
row.names(col_annos) <- col_annos$geneSymbol


## get speckle protein gene expression
expressionSpeckleSignature <- expression[col_annos$geneSymbol,2:312]
expressionSpeckleSignature <- na.omit(expressionSpeckleSignature)
expressionSpeckleSignature <- t(expressionSpeckleSignature) # so that rows are samples, columns are gene
expressionSpeckleSignature <- as.data.frame(expressionSpeckleSignature)

## get z scores
N <- scale(expressionSpeckleSignature, center=TRUE, scale=TRUE)
N <- as.data.frame(N)

## calculate speckle score, which is the sum of the z scores of Igenes * Imultiplier + sum of the z scores of IIgenes * IImultiplier
N$speckleScore <- rowSums(N[,col_annos$geneSymbol[col_annos$signature=="I"]]*Imultiplier) + rowSums(N[,col_annos$geneSymbol[col_annos$signature=="II"]]*IImultiplier)
N$arm <- metadata$Arm[match(row.names(N), metadata$RNA_ID)]
N$benefit <- metadata$Benefit[match(row.names(N), metadata$RNA_ID)]
N$tumorShrinkage <- metadata$Tumor_Shrinkage[match(row.names(N), metadata$RNA_ID)]
N$MSKCC <- metadata$MSKCC[match(row.names(N), metadata$RNA_ID)]
N$PFS <- metadata$PFS[match(row.names(N), metadata$RNA_ID)]
N$PFS_CNSR <- metadata$PFS_CNSR[match(row.names(N), metadata$RNA_ID)]
N$OS <- metadata$OS[match(row.names(N), metadata$RNA_ID)]
N$OS_CNSR <- metadata$OS_CNSR[match(row.names(N), metadata$RNA_ID)]

N$signature <- "NA"
N$signature[N$speckleScore > 0] <- "I"
N$signature[N$speckleScore < 0] <- "II"
N$signatureNtile <- ntile(N$speckleScore, 3)
N$signatureNtile <- as.factor(N$signatureNtile)

## checking different variations of survival curves ##
# survival by signature, combined treatments
sfit <- survfit(Surv(OS, OS_CNSR)~signature, data = N)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))
# survival by signature nivolumab
sfit <- survfit(Surv(OS, OS_CNSR)~signature, data = N[N$arm == "NIVOLUMAB",])
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))
# survival by signature everolimus
sfit <- survfit(Surv(OS, OS_CNSR)~signature, data = N[N$arm == "EVEROLIMUS",])
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))

# progression-free survival by signature and treatment
sfit <- survfit(Surv(PFS, PFS_CNSR)~arm, data = N[N$signature == "I",])
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))
sfit <- survfit(Surv(PFS, PFS_CNSR)~arm, data = N[N$signature == "II",])
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))

# overall survival by signature and treatment with graphs saved, and sample sizes printed
sfit <- survfit(Surv(OS, OS_CNSR)~arm, data = N[N$signature == "I",])
p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))
filename = "OS_signatureI.pdf"
pdf(filename, width = 4.5, height = 5, onefile=FALSE)
print(p)
dev.off()
length(row.names(N[N$signature == "I" & N$arm == "NIVOLUMAB",]))
length(row.names(N[N$signature == "I" & N$arm == "EVEROLIMUS",]))
sfit <- survfit(Surv(OS, OS_CNSR)~arm, data = N[N$signature == "II",])
p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE,  palette = c("aquamarine4", "grey1"), xlim = c(0,60))
filename = "OS_signatureII.pdf"
pdf(filename, width = 4.5, height = 5, onefile=FALSE)
print(p)
dev.off()
length(row.names(N[N$signature == "II" & N$arm == "NIVOLUMAB",]))
length(row.names(N[N$signature == "II" & N$arm == "EVEROLIMUS",]))
