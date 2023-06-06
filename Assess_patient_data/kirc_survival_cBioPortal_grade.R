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
patient.data <- patient.data[patient.data$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE I" | patient.data$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE II",]

## add binary for speckle signature
patient.data$SpeckleSignature <- with(patient.data, ifelse(SpeckleScore > 0, "I", "II"))


if (!dir.exists("cBioPortal_survival/survival_byGrade")){
  dir.create("cBioPortal_survival/survival_byGrade")
}

## plot survival based on Sig I or II
dssfit <- survfit(Surv(DSS_MONTHS, as.numeric(str_sub(patient.data$DSS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = patient.data)
p = ggsurvplot(dssfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I", "II"), legend.title="Speckle signature", palette = c("#c66b3d", "#26495c"))
fileName = paste("cBioPortal_survival/survival_byGrade/DSS_", study, "_IvsII_stageIandII.pdf", sep = "")
pdf(fileName, width = 5, height = 5, onefile=FALSE)
print(p)
dev.off()


## same as above, but with stage III and IV
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

patient.data <- patient.data[patient.data$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE III" | patient.data$AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE IV",]

## add binary for speckle signature
patient.data$SpeckleSignature <- with(patient.data, ifelse(SpeckleScore > 0, "I", "II"))

dssfit <- survfit(Surv(DSS_MONTHS, as.numeric(str_sub(patient.data$DSS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = patient.data)
p = ggsurvplot(dssfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I", "II"), legend.title="Speckle signature", palette = c("#c66b3d", "#26495c"))
fileName = paste("cBioPortal_survival/survival_byGrade/DSS_", study, "_IvsII_stageIIIandIV.pdf", sep = "")
pdf(fileName, width = 5, height = 5, onefile=FALSE)
print(p)
dev.off()



