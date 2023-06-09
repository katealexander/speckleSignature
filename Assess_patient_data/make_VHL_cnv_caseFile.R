setwd("~/Documents/speckleSignature/Assess_patient_data/")

library(stringr)

# change this value to the path of the directory of the CNV data
cnv.data.dir <- "~/Desktop/speckleSignatureAnalysis/gdc_kirc_CNV_20230609/"

# load sample sheet, which links samples to file names
sample.sheet <- read.table("gdc_sample_sheet_CNV.2023-06-08.tsv", sep = "\t", header = T)


case.VHL.cnv <- data.frame(id = character(0), VHLminCopyNumber = numeric(0))

idsLength <- length(list.dirs(path = cnv.data.dir, full.names = F))

for (id in list.dirs(path = cnv.data.dir, full.names = F)[2:idsLength]){
  id
  fname <- list.files(paste(cnv.data.dir, id, sep = ""))[grep(pattern = "TCGA", list.files(paste(cnv.data.dir, id, sep = "")))]
  full.path <- paste(cnv.data.dir, id, "/", fname, sep="")
  cnv.data <- read.table(full.path, sep = "\t", header=T)
  min.copy.number.VHL <- cnv.data[cnv.data$gene_name == "VHL",]$min_copy_number
  caseID <- strsplit(sample.sheet$Case.ID[sample.sheet$File.Name == fname], ",")[[1]][1]
  case.VHL.cnv <- rbind(case.VHL.cnv, setNames(as.list(c(caseID, min.copy.number.VHL)), names(case.VHL.cnv)))
}
write.table(case.VHL.cnv, row.names = F, quote = F, file = "KIRC_VHL_CNV.txt", sep = "\t")

