setwd("~/Documents/speckleSignature/Assess_patient_data/")

library(stringr)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(hrbrthemes)
library(viridis)

study <- "kirc_tcga_pan_can_atlas_2018"

## get clinical data
patient.data <- read.table(paste("patientData_cBioPortal/", study, "_patientData.txt", sep=""), sep="\t", header=T)

## get speckle scores
speckle.data <- read.table(paste("../Speckle_signature_definition/cBioPortal_signatureSpeckleScores_22cancers/", study, "_Zscore_speckleScore.txt", sep=""), sep="\t", header=T, row.names=1)
speckle.data.scores <- as.data.frame(speckle.data$speckleScore)
row.names(speckle.data.scores) <- row.names(speckle.data)
colnames(speckle.data.scores) <- "SpeckleScore"
speckle.data.scores$patientId <- str_sub(row.names(speckle.data.scores), start = 1, end = 12)
speckle.data.scores$patientId <- gsub("\\.", "-", speckle.data.scores$patientId)

## add speckle scores to clinical data
patient.data <- merge(patient.data, speckle.data.scores, by = "patientId")

## add binary for speckle signature
patient.data$SpeckleSignature <- with(patient.data, ifelse(SpeckleScore > 0, "I", "II"))

## get VHL point mutant samples or samples with CNV loss
VHLmut.samples <- read.table("KIRC_TCGA_VHL_casesEffected.2023-02-24.tsv", header = T)
VHLcnv.status <- read.table("KIRC_VHL_CNV.txt", header = T)
VHLloss.samples <- as.data.frame(VHLcnv.status$id[VHLcnv.status$VHLminCopyNumber == 1])
VHLnorm.samples <- as.data.frame(VHLcnv.status$id[VHLcnv.status$VHLminCopyNumber == 2])
colnames(VHLloss.samples) <- "id"
colnames(VHLnorm.samples) <- "id"
patient.data.VHLmut <- patient.data[toupper(patient.data$OTHER_PATIENT_ID) %in% intersect(toupper(VHLmut.samples$id), toupper(patient.data$OTHER_PATIENT_ID)),]
patient.data.VHLcnvLoss <- patient.data[toupper(patient.data$patientId) %in% intersect(toupper(VHLloss.samples$id), toupper(patient.data$patientId)),]
patient.data.VHLcnvNorm <- patient.data[toupper(patient.data$patientId) %in% intersect(toupper(VHLnorm.samples$id), toupper(patient.data$patientId)),]
patient.data.VHL.LOF <- unique(rbind(patient.data.VHLmut, patient.data.VHLcnvLoss))

patient.data.VHLwt <- patient.data.VHLcnvNorm[toupper(patient.data.VHLcnvNorm$OTHER_PATIENT_ID) %in% setdiff(toupper(patient.data.VHLcnvNorm$OTHER_PATIENT_ID), toupper(patient.data.VHL.LOF$OTHER_PATIENT_ID)),]

if (!dir.exists("cBioPortal_survival/survival_KIRC_VHLmutStatus")){
  dir.create("cBioPortal_survival/survival_KIRC_VHLmutStatus")
}

## plot survival based on Sig I or II
dssfit <- survfit(Surv(DSS_MONTHS, as.numeric(str_sub(patient.data.VHLmut$DSS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = patient.data.VHLmut)
p = ggsurvplot(dssfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I", "II"), legend.title="Speckle signature", palette = c("#c66b3d", "#26495c"))
fileName = paste("cBioPortal_survival/survival_KIRC_VHLmutStatus/DSS_", study, "_IvsII_KIRCmut.pdf", sep = "")
pdf(fileName, width = 5, height = 5, onefile=FALSE)
print(p)
dev.off()

dssfit <- survfit(Surv(DSS_MONTHS, as.numeric(str_sub(patient.data.VHLwt$DSS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = patient.data.VHLwt)
p = ggsurvplot(dssfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I", "II"), legend.title="Speckle signature", palette = c("#c66b3d", "#26495c"))
fileName = paste("cBioPortal_survival/survival_KIRC_VHLmutStatus/DSS_", study, "_IvsII_KIRCwt.pdf", sep = "")
pdf(fileName, width = 5, height = 5, onefile=FALSE)
print(p)
dev.off()



