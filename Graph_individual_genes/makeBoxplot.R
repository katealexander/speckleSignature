setwd("~/Documents/speckleSignature/Graph_individual_genes/")
library(ggplot2)
library(ggpubr)


for (file in list.files(path = "individualGenes/")){
  if (grepl("_patientGroups.txt", file)){
    A <- read.table(paste("individualGenes/", file, sep = ""), sep="\t")
    gene <- strsplit(file, split = "_patientGroups")[[1]][1]
    # my_specific_comparisons <- list(c("group1", "group2"), c("group2", "group3"), c("group1", "group3"))
    my_specific_comparisons <- list(c("group1", "group2"), c("group1", "group3"))
    p = ggplot(A, aes(x=V2, y=V1, fill=factor(V2))) +
      geom_boxplot(col="black", size=1, outlier.size = 0) +
      geom_jitter(color="grey20", size=.25, alpha= .5) +
      theme_classic() +
      theme(legend.position="none") +
      stat_compare_means(comparisons = my_specific_comparisons, method = "wilcox.test", size = 4, tip.length = 0.015, step.increase = 0.08, label = "p.signif", vjust = .6, hide.ns = T, ) +
      scale_x_discrete(labels = c("I", "N", "II")) + 
      xlab("KIRC Patient Group") +
      ylab(paste(gene, "normalized counts")) +
      ylim(0,NA) +
      scale_fill_manual(values = c("#D4A78E", "#BCBEC0", "#9CACB4"))
    fileName = paste("individualGenes/", gene, "_boxplot.pdf", sep="")
    pdf(fileName, width = 2, height = 2.25, onefile=FALSE)
    print(p)
    dev.off()
  }
}

