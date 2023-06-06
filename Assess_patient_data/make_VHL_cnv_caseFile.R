setwd("~/Documents/speckleSignature/Assess_patient_data/")

library(stringr)


case.VHL.cnv <- data.frame(id = character(0), VHLminCopyNumber = numeric(0))

for (id in list.dirs(path = "gdc_kirc_CNV", full.names = F)[2:length(list.dirs(path = "gdc_kirc_CNV", full.names = F))]){
  if (length(list.files(paste("gdc_kirc_CNV/", id, sep=""))) == 2) {
    annotation <- read.table(paste("gdc_kirc_CNV/", id, "/annotations.txt", sep = ""), header=T, sep = "\t")
    caseID <- annotation$entity_id[1]
    file.name <- list.files(paste("gdc_kirc_CNV/", id, sep=""))[grep(pattern = "TCGA", list.files(paste("gdc_kirc_CNV/", id, sep="")))]
    full.path <- paste(getwd(), "/gdc_kirc_CNV/", id, "/", file.name, sep = "")
    cnv.data <- read.table(full.path, sep = "\t", header=T)
    min.copy.number.VHL <- cnv.data[cnv.data$gene_name == "VHL",]$min_copy_number
    case.VHL.cnv <- rbind(case.VHL.cnv, setNames(as.list(c(caseID, min.copy.number.VHL)), names(case.VHL.cnv)))
  }
}
write.table(case.VHL.cnv, row.names = F, quote = F, file = "KIRC_VHL_CNV.txt", sep = "\t")

