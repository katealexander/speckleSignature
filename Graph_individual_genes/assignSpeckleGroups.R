setwd("~/Documents/speckleSignature/Graph_individual_genes/")

library(stringr)
library(dplyr)

study <- "kirc_tcga_pan_can_atlas_2018"

## get speckle scores
speckle.data <- read.table(paste("../Speckle_signature_definition/cBioPortal_signatureSpeckleScores_22cancers/", study, "_Zscore_speckleScore.txt", sep=""), sep="\t", header=T, row.names=1)
speckle.data.scores <- as.data.frame(speckle.data$speckleScore)
row.names(speckle.data.scores) <- row.names(speckle.data)
colnames(speckle.data.scores) <- "SpeckleScore"
speckle.data.scores$patientId <- str_sub(row.names(speckle.data.scores), start = 1, end = 12)
speckle.data.scores$patientId <- gsub("\\.", "-", speckle.data.scores$patientId)
speckle.data.scores$ntile <- ntile(speckle.data.scores$SpeckleScore, 4)

## get VHL mutant samples
sample.sheet <- read.table("gdc_sample_sheet_KIRC_VHLmut.2021-09-20.tsv", sep = "\t", header = T)

## assign the speckle groups 1 is signature I, 3 is Signature II, 2 is normal
sample.sheet$speckleGroup <- NA
sample.sheet$speckleGroup[sample.sheet$Case.ID %in% speckle.data.scores$patientId[speckle.data.scores$ntile == "1"]] <- "3"
sample.sheet$speckleGroup[sample.sheet$Case.ID %in% speckle.data.scores$patientId[speckle.data.scores$ntile == "4"]] <- "1"
sample.sheet$speckleGroup[sample.sheet$Sample.Type == "Solid Tissue Normal"] <- "2"
sample.sheet <- na.omit(sample.sheet)

## write the file
write.table(as.data.frame(sample.sheet[, c(7,9)]), file="KIRC_VHLmutant_specklePatientGroups.txt", sep = "\t", quote = F, row.names = F)
 
