# Graph individual genes
This pipeline provides instructions for making boxplots in R graphing individual gene expression from the KIRC TCGA cohort split into samples that were tumor speckle Signature I, tumor speckle Signature II, or tumor-adjacent "normal" tissue. 

# Gene expression file
RNA expression STAR gene counts TCGA-KIRC project were downloaded from [GDC](https://portal.gdc.cancer.gov/), formatted into table format with FPKM-UQ normalized expression values:

```python makeExpressionTable_GDC.py gdc_sample_sheet.RNA.2023-06-09.tsv ENS_geneName.txt ~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/gdc_kirc_RNA_20230609_191255.883495/*/*.tsv > ~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/KIRC_FPKM_UQ.txt```

The expression files are not included in the repository. The resulting expression file contains ENSEMBL gene id, Gene Symbol, then sampleIDs as columns and genes as rows. 

# Assign speckle groups
The following script calculates speckle scores and assigns samples to Signature I (Group 1), normal adjacent (Group 2), or Signature II (Group 3). Here, Signature I and Signature II are defined as the top and bottom quarter of speckle signature scores, respectively. 

This script has the option to filter samples based on VHL expression status. See comments within script.

```Rscript assignSpeckleGroups.R```

The above script will output the file, "KIRC_specklePatientGroups.txt".

# Format gene expression for plotting
The following script will extract data from an individual gene and format it for making a boxplot in R. In this example, a list of genes, "genesOfInterest.txt" is used. 

```while IFS= read -r line; do python getExpressionForPatientGroupsOfGenes.py KIRC_specklePatientGroups.txt ~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/KIRC_FPKM_UQ.txt $line > "individualGenes/"$line"_patientGroups.txt"; done < genesOfInterest.txt```

The above script outputs formated gene-specific text files into the folder, "individualGenes"

# Plot gene expression by speckle group
To generate boxplots for the above genes:

```Rscript makeBoxplot.R```

The above script generates boxplots where each dot represents a sample. I - Signature I, N - normal adjacent tissue, II - Signature II

# Generate table of median speckle patient group expression per gene
```python getMedianExpressionOfPatientGroup.py hg19_genes.txt KIRC_specklePatientGroups.txt ~/Desktop/Alexander2023_filesTooBigForGithub/speckleSignature/KIRC_FPKM_UQ.txt > medianGeneExpression_KIRC_specklepatientGroups.txt```

The above script also performs Student's t-tests (not corrected) for differences between each Signature group and normal and between Signature I and II.

# Assess HIF2A target genes in the KIRC patient groups
The following script returns a heatmap of the median HIF2A target gene expression within the patient speckle groups and normal adjacent tissue. It also performs cluster compare to assess the functional categories of Signature I-high HIF2A target genes versus Signature II-high and non differential HIF2A targets. 

Here, HIF2A target genes are defined by the following: 
1. Decreasing upon treatment with a HIF2A inhibitor, PT2399, in 786O ccRCC cell line or Increasing upon induced HIF2A expression in MCF7 cells (gene list in file "HIF2Atargets_MCF7_786O_combined.txt")
2. AND Higher expression in ccRCC tumor compared to normal in TCGA data. 

```Rscript HIFresponsiveGenePlots.R```

The above script returns two files:
1. medianExpression_specklePatientGroups_HIFresponsiveGenes_heatmap.pdf - a heatmap of the HIF2A target genes' median expression in Signature I, Normal, and Signature II groups
2. clusterProfiler_IvsIIbias_withNS.pdf - functional analysis of I-biased, non-differential, and II-biased HIF2A target genes.












