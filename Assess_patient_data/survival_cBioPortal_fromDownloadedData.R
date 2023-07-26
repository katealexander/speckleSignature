
setwd("~/Documents/speckleSignature/Assess_patient_data/")

library(stringr)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(stats)



## create new directories to store files if they don't already exist
if (!dir.exists("cBioPortal_survival")){
  dir.create("cBioPortal_survival")
}

if (!dir.exists("cBioPortal_survival/survival_IvsII")){
  dir.create("cBioPortal_survival/survival_IvsII")
}

if (!dir.exists("cBioPortal_survival/survival_topQuarterScores")){
  dir.create("cBioPortal_survival/survival_topQuarterScores")
}

# if (!dir.exists("cBioPortal_survival/continuousData_speckleNtile_boxplots")){
#   dir.create("cBioPortal_survival/continuousData_speckleNtile_boxplots")
# }

if (!dir.exists("cBioPortal_survival/categoricalData_speckleNtile_bargraphs")){
  dir.create("cBioPortal_survival/categoricalData_speckleNtile_bargraphs")
}



for (file in list.files(path = "./patientData_cBioPortal/")){
  if (grepl("_patientData.txt", file)){
    study <- strsplit(file, split = "_patientData.txt")[[1]][1]
    # study <- "brca_tcga_pan_can_atlas_2018" # for testing, should be commented out
    
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
    
    ## Add column for speckle ntiles
    patient.data$SpeckleNtile <- as.factor(ntile(patient.data$SpeckleScore, 4))
    
    ## plot survival based on Sig I or II
    if ("DSS_MONTHS" %in% colnames(patient.data)){
      dssfit <- survfit(Surv(DSS_MONTHS, as.numeric(str_sub(patient.data$DSS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = patient.data)
      p = ggsurvplot(dssfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I", "II"), legend.title="Speckle signature", palette = c("#c66b3d", "#26495c"))
      fileName = paste("cBioPortal_survival/survival_IvsII/DSS_", study, "_IvsII.pdf", sep = "")
      pdf(fileName, width = 5, height = 5, onefile=FALSE)
      print(p)
      dev.off()
    }
    
    osfit <- survfit(Surv(OS_MONTHS, as.numeric(str_sub(patient.data$OS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = patient.data)
    p = ggsurvplot(osfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I", "II"), legend.title="Speckle signature", palette = c("#c66b3d", "#26495c"))
    fileName = paste("cBioPortal_survival/survival_IvsII/OS_", study, "_IvsII.pdf", sep = "")
    pdf(fileName, width = 5, height = 5, onefile=FALSE)
    print(p)
    dev.off()

    
    ## plot survival based on top and bottom 25%
    extremes <- patient.data[(patient.data$SpeckleNtile == "1") | (patient.data$SpeckleNtile == "4"),]
    if ("DSS_MONTHS" %in% colnames(patient.data)){
      dssfit <- survfit(Surv(DSS_MONTHS, as.numeric(str_sub(extremes$DSS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = extremes)
      p = ggsurvplot(dssfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I top 25%", "II top 25%"), legend.title="Speckle Score", palette = c("#c66b3d", "#26495c"))
      fileName = paste("cBioPortal_survival/survival_topQuarterScores/DSS_", study, "_topQuarterScores.pdf", sep = "")
      pdf(fileName, width = 5, height = 5, onefile=FALSE)
      print(p)
      dev.off()
    }
    
    osfit <- survfit(Surv(OS_MONTHS, as.numeric(str_sub(extremes$OS_STATUS, start = 1, end = 1)))~SpeckleSignature, data = extremes)
    p = ggsurvplot(osfit, conf.int=TRUE, pval=TRUE, legend.labs=c("I top 25%", "II top 25%"), legend.title="Speckle Score", palette = c("#c66b3d", "#26495c"))
    fileName = paste("cBioPortal_survival/survival_topQuarterScores/OS_", study, "_topQuarterScores.pdf", sep = "")
    pdf(fileName, width = 5, height = 5, onefile=FALSE)
    print(p)
    dev.off()

    
    # ## Boxplots of continuous data for each speckle decile
    # continuous.data.headers <- c("BUFFA_HYPOXIA_SCORE", "RAGNUM_HYPOXIA_SCORE", "WINTER_HYPOXIA_SCORE", "ANEUPLOIDY_SCORE", "FRACTION_GENOME_ALTERED", "MSI_SCORE_MANTIS",
    #                              "MSI_SENSOR_SCORE", "MUTATION_COUNT", "AGE", "TMB_NONSYNONYMOUS")
    # mycomparisons <- list(c("1","4"))
    # for (contData in continuous.data.headers){
    #   if (contData %in% colnames(patient.data)){
    #     ## plot the continuous data
    #     p = ggplot(patient.data, aes(x=SpeckleNtile, y=patient.data[[contData]], fill=factor(SpeckleNtile))) + 
    #       geom_boxplot(col="black", size=1, outlier.size = 0) +
    #       geom_jitter(color="grey20", size=.5, alpha= .5) +
    #       theme_classic() +
    #       scale_fill_viridis(discrete = TRUE, alpha=0.8, option = "B") +
    #       theme(legend.position="none") +
    #       stat_compare_means(comparisons = mycomparisons, method = "wilcox.test", size = 4, tip.length = 0.02) +
    #       ylab(contData)
    #     fileName = paste("cBioPortal_survival/continuousData_speckleNtile_boxplots/", contData, "_", study, ".pdf", sep = "")
    #     pdf(fileName, width = 3, height = 4, onefile=FALSE)
    #     print(p)
    #     dev.off()
    #   }
    # }
    
  
    ## Stacked bargraphs of categorical data for each speckle decile
    categorical.data.headers <- c("AJCC_PATHOLOGIC_TUMOR_STAGE", "ETHNICITY", "HISTORY_NEOADJUVANT_TRTYN", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
                                  "PATH_M_STAGE", "PATH_N_STAGE", "PATH_T_STAGE", "PERSON_NEOPLASM_CANCER_STATUS", "PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT",
                                  "PRIOR_DX", "RACE", "RADIATION_THERAPY", "SEX", "GRADE", "TISSUE_PROSPECTIVE_COLLECTION_INDICATOR", "TISSUE_SOURCE_SITE_CODE", "SUBTYPE",
                                  "CANCER_TYPE_DETAILED")
    for (catData in categorical.data.headers){
        # catData <- "RACE"
        if (catData %in% colnames(patient.data)){
          x <- table(patient.data[[catData]], patient.data$SpeckleNtile)
          x <- as.matrix(x)
          if (length(row.names(x)) >1){
            p = fisher.test(x, simulate.p.value = TRUE, B = 500000)$p.value
            p <- formatC(p, format="e", digits = 2)
            x <- as.data.frame(x)
            colnames(x) <- c(catData,"SpeckleScoreNtile","Frequency")
            p = ggplot(x, aes(fill=x[[catData]], y=Frequency, x=SpeckleScoreNtile)) +
              geom_bar(position="stack", stat="identity") +
              scale_fill_viridis(discrete = T) +
              theme_classic() +
              xlab("Speckle Score Ntile") +
              labs(fill = "Category", title = catData, subtitle = paste("Fishers Exact p-value: ", p, sep = ""))
            fileName = paste("cBioPortal_survival/categoricalData_speckleNtile_bargraphs/", catData, "_", study, ".pdf", sep = "")
            pdf(fileName, width = 4, height = 4, onefile=FALSE)
            print(p)
            dev.off()
        }
      }
    } 
    
  }
  
}
  