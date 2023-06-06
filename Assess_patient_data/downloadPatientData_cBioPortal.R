## This script accesses TCGA data via cBioPortal

setwd("~/Documents/speckleSignature/Assess_patient_data/")

library(cBioPortalData)
library(AnVIL)



##### Access to cBioPortal API for R #####

## cBioPortal TCGA API is used to generate the speckle signature from tumor samples
cbio <- cBioPortal()

## Table of studies with RNA previously generated in dataAccess_cBioPortal.R
studiesWithRNA <- read.table("../Speckle_signature_definition/studiesWithRNA_cBioPortal.txt", sep = "\t", header = T)


##### Loop through the cBioPortal TCGA studies and get patient data #####

## make folders to hold files created
if (!dir.exists("patientData_cBioPortal")){
  dir.create("patientData_cBioPortal")
}

for (id in studiesWithRNA$studyId){
  patient.data <- clinicalData(cbio, id)
  write.table(patient.data, sep="\t", file=paste("patientData_cBioPortal/", id, "_patientData.txt", sep=""), quote = F, row.names = F)
}


