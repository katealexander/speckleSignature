## This script accesses TCGA data via cBioPortal

setwd("~/Documents/speckleSignature/Speckle_signature_definition/")

library(cBioPortalData)
library(AnVIL)



##### Access to cBioPortal API for R #####

## cBioPortal TCGA API is used to generate the speckle signature from tumor samples
cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)
rna.experiment.name <- "mrna_seq_v2_rsem"



##### Getting speckle protein genes #####

## subcellular localization data downloaded from Human Protein atlas on 02/12/2023
subcellular <- read.table("subcellular_location.tsv", sep = "\t", header =T)
s = "Nuclear speckles"

## get the proteins that are annotated as Enhanced, Supported, or Approved in nuclear speckles
subcellularSpeckle <- subcellular[grep(s, subcellular$Enhanced),]
subcellularSpeckle <- rbind(subcellularSpeckle, subcellular[grep(s, subcellular$Supported),])
subcellularSpeckle <- rbind(subcellularSpeckle, subcellular[grep(s, subcellular$Approved),])
row.names(subcellularSpeckle) <- subcellularSpeckle$Gene.name


#### Defining the studies to access #####

n <- 50 # minimum number of cases in study
studiesToExclude <- c("pcpg_tcga_pub", "stad_tcga_pub", "blca_tcga_pub_2017", "blca_tcga_pub", "brca_tcga_pub2015", "gbm_tcga_pub2013", "hnsc_tcga_pub", "kich_tcga_pub", "kirc_tcga_pub", "laml_tcga_pub", "luad_tcga_pub", "pcpg_tcga_pub", "prad_tcga_pub", "sarc_tcga_pub", "thca_tcga_pub", "ucec_tcga_pub") ## excluded studies because of redundancy with pancancer datasets
studiesToAssess <- studies
studiesToAssess <- studies[!(studies$studyId %in% studiesToExclude),]
studiesToAssess <- studiesToAssess[studiesToAssess$api_build == "TRUE",] # has api build available
studiesToAssess <- na.omit(studiesToAssess)
studiesToAssess <- studiesToAssess[studiesToAssess$allSampleCount > n,]
studiesToAssess <- studiesToAssess[grep("tcga", studiesToAssess$studyId),]
#studiesToAssess <- studiesToAssess[35:51,] # for testing, comment out to look at all studies


##### Loop through the cBioPortal TCGA studies and get rotations for speckle protein genes #####

## dataframe to keep track of which studies have RNA data
studiesWithRNA <- data.frame(matrix(ncol = 15, nrow = 0)) 
colnames(studiesWithRNA) <- colnames(studiesToAssess)

## make folders to hold files created
if (!dir.exists("speckleProteinGeneExpression")){
  dir.create("speckleProteinGeneExpression")
}

for (id in studiesToAssess$studyId){
  ## check for expression data
  mols <- molecularProfiles(cbio, id)
  mols <- as.data.frame(mols)
  RNAprofileID <- mols$molecularProfileId[grep("_rna_seq_v2_mrna$", mols$molecularProfileId)]
  
  if (!identical(RNAprofileID, character(0))){
    studiesWithRNA <- rbind(studiesWithRNA, studiesToAssess[studiesToAssess$studyId == id,]) # adds study to the list of studies with expression data
    expressionTibble <- getDataByGenes(cbio, molecularProfileIds = RNAprofileID, sampleIds = allSamples(cbio, id)$sampleId, genes = subcellularSpeckle$Gene.name, by = 'hugoGeneSymbol')
    expressionTibble <- na.omit(expressionTibble)
    expressionDF <- data.frame(matrix(ncol = length(unique(expressionTibble[[RNAprofileID]]$sampleId)), nrow = length(unique(expressionTibble[[RNAprofileID]]$hugoGeneSymbol))))
    colnames(expressionDF) <- unique(expressionTibble[[RNAprofileID]]$sampleId)
    row.names(expressionDF) <- unique(expressionTibble[[RNAprofileID]]$hugoGeneSymbol)
    
    for (i in colnames(expressionDF)){
      for (j in row.names(expressionDF)){
        value = expressionTibble[[RNAprofileID]]$value[expressionTibble[[RNAprofileID]]$sampleId == i & expressionTibble[[RNAprofileID]]$hugoGeneSymbol == j]
        if (!identical(value, numeric(0))){
          expressionDF[j, i] <- as.numeric(value)
        }
      }
    }
    
    expressionDF <- expressionDF[rowMedians(as.matrix(expressionDF)) > 150,]
    speckle.expression.data <- na.omit(expressionDF)
    write.table(speckle.expression.data, sep="\t", file=paste("speckleProteinGeneExpression/", id, "_speckleExpression.txt", sep=""), quote = F)

    ## pca on speckle protein gene expression, and extract rotations
    pca <- prcomp(t(speckle.expression.data), center=TRUE, scale. = TRUE)
    rotations <- as.data.frame(pca$rotation)
    PC1 <- as.data.frame(rotations$PC1)
    row.names(PC1) <- row.names(rotations)
    colnames(PC1) <- "PC1"
    
    ## if SON has a positive rotation, reverse the signs so that cancers can be compared
    ## Sig I will be positive; Sig II, which includes SON will be negative
    SON.PC1 <- rotations["SON", "PC1"] # get the SON PC1 rotation
    if (SON.PC1 > 0) {
      PC1$PC1 <- -PC1$PC1
      }
    
    ## add study rotations to the subcellularSpeckle dataframe
    PC1$Gene.name <- row.names(PC1)
    colnames(PC1) <- c(id, "Gene.name")
    subcellularSpeckle <- merge(subcellularSpeckle, PC1, all=TRUE, by="Gene.name")
  }
}

write.table(studiesWithRNA, sep="\t", file="studiesWithRNA_cBioPortal.txt", quote = F)
write.table(subcellularSpeckle, sep="\t", file="speckleProteinGenes_rotations_allStudiesWithRNA.txt", quote = F)

######### END ##########
