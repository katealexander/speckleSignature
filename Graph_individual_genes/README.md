# Graph individual genes
This pipeline provides instructions for making boxplots in R graphing individual gene expression from the KIRC VHL mutant TCGA cohort split into samples that were tumor speckle Signature I, tumor speckle Signature II, or tumor-adjacent "normal" tissue. 

# GDC downloads
The following was downloaded from the [GDC](https://portal.gdc.cancer.gov/) on 09/20/2021:
1. "gdc_sample_sheet_KIRC_VHLmut.2021-09-20.tsv" - The VHL-mutant KIRC TCGA cohort sample sheet containing gene expression File ID and corresponding Case ID, Sample ID, and Sample Type. Included in this repository.
2. RNA expression data (FPKM-UQ) files of the samples from "gdc_sample_sheet_KIRC_VHLmut.2021-09-20.tsv". Not included in this repository.

# Formatting expression data
I used a Python script to format the expression data downloaded from the GDC into a table format where rows were genes, the first column was Ensembl gene names, the second column was Official Gene Symbols, and the remaining columns were TCGA samples.

```python makeExpressionTable_GDC.py gdc_sample_sheet_KIRC_VHLmut.2021-09-20.tsv ENS_geneName.txt pathToExpression/*/*UQ.txt > KIRC_VHLmutant_FPKM_UQ.txt```

The file, "ENS_geneName.txt", included in this repository, is a file linking Ensembl gene names to gene symbols downloaded from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) on 09/21/2021

# Assign speckle groups
The following script uses the speckle signature scores generated in "Speckle_signature_definition" and assigns the VHL mutant samples to Signature I (Group 1), normal adjacent (Group 2), or Signature II (Group 3). Here, Signature I and Signature II are defined as the top and bottom 25% of speckle signature scores, respectively.

```Rscript assignSpeckleGroups.R```

The above script will output the file, "KIRC_VHLmutant_specklePatientGroups.txt".

# Format gene expression for plotting
The following script will extract data from an individual gene and format it for making a boxplot in R. In this example, a list of genes, "genesOfInterest.txt" is used. 

```while IFS= read -r line; do python getExpressionForPatientGroupsOfGenes.py KIRC_VHLmutant_specklePatientGroups.txt KIRC_VHLmutant_FPKM_UQ.txt $line > "individualGenes/"$line"_patientGroups.txt"; done < genesOfInterest.txt```

The above script outputs formated gene-specific text files into the folder, "individualGenes"

# Plot gene expression by speckle group
To generate boxplots for the above genes:

```Rscript makeBoxplot.R```

The above script generates boxplots where each dot represents a sample. I - Signature I, N - normal adjacent tissue, II - Signature II

# Generate table of median speckle patient group expression per gene
```python getMedianExpressionOfPatientGroup.py hg19_genes.txt KIRC_VHLmutant_specklePatientGroups.txt KIRC_VHLmutant_FPKM_UQ.txt > medianGeneExpression_KIRC_VHLmutant_specklepatientGroups.txt```

The above script also performs Student's t-tests (not corrected) for differences between each Signature group and normal and between Signature I and II.









