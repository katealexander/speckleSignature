# Defining the speckle signature
This pipeline describes how I arrived at the speckle protein genes that comprise the RNA-seq based speckle signature

# Speckle-resident proteins
The analysis below is based on speckle-residency annotations the file "subcellular_location.tsv", downloaded from the [The Human Protein Atlas](https://www.proteinatlas.org/about/download) on 02/12/2023. Speckle-resident proteins were defined as those annotated as "Enhanced", "Supported", or "Approved" for subcellular localization in nuclear speckles (446 speckle proteins). 

# Getting PC1 rotations of speckle protein genes
### Rationale
My goal was to determine how speckle protein gene expression varied in cancer. To extract which speckle protein genes contribute the most to patient variation, I performed Principal Comonent Analysis using the prcomp function in R with speckle protein gene expression as input. In Principal Component Analysis, each speckle protein gene is assigned a rotation that separates samples along each Principal Component, with Prinicpal Component 1 (PC1) reflecting the most patient variation. The absolute value and direction (positive or negative) of speckle protein gene rotations thus reflect the contributions of that speckle protein gene to patient variation.
### Making PC1 rotation signs consistent between cancer types
In previous analyses, I noticed that the top-contributing speckle protein genes (high absolute value PC1 rotation) fall into two reciprocally-expressed groups of speckle protein genes that have opposite PC1 rotation signs (positive versus negative). The speckle marker protein, SON, consistently resided within one of these groups. To make the PC1 rotation sign consistent across cancer types, I required the SON PC1 rotation to be negative and flipped PC1 rotation signs if it was not. Speckle protein genes with positive PC1 rotations are high in speckle Signature I; speckle protein genes with negative PC1 rotations are high in speckle Signature II.
### Gene expression cutoff
For a more robust analysis, I included only speckle protein genes with median expression of >150 rsem. This was assessed separately for each cancer type.
### Sample size cutoff
I excluded studies with fewer than 50 samples.
### Getting a table of speckle protein gene rotations
The following R script extracts expression of speckle protein genes from cBioPortal
```Rscript dataAccess_cBioPortal.R```
It returns:
1. "speckleProteinGeneExpression" - A folder containing expression of speckle protien gene for each cancer study downloaded
2. "speckleProteinGenes_rotations_allStudiesWithRNA.txt" - A table of PC1 rotations for speckle protein genes for each study. Speckle protein genes with median expression lower than rsem of 150 will be "NA"
3. "studiesWithRNA_cBioPortal.txt", which is the subset of datsets with RNA expression data available (mrna_seq_v2_rsem) for download via cBioPortal
Each of these are included in this repository, downloaded 06/05/2023. 

# Comparing speckle protein gene PC1 rotations between cancer types
The following R script takes the rotations generated in "dataAccess_cBioPortal.R", calculates Pearson correlations, and returns a heatmap of the pairwise Pearson correlations
```Rscript rotationCorrelations.R```
It shows that 24 of 30 cancer types have highly correlated speckle protein gene PC1 rotations, with 22 of these having particularly strong correlations. This indicates consistent speckle protein gene contributions to patient variation across multiple cancer types, justifying the use of a speckle signature definition that can be used across cancer types.

# Extracting consistent speckle protein genes for speckle signature
Over the years, I have tried a few different methods for deciding which speckle protein genes to use for defining the multi-cancer speckle signature. Each method has given similar results, indicating a robustness to variation in the ultimate selection of speckle protein genes. The below method counts up how many studies of the 22 highly-consistent cancer types above have positively (Signature I high) or negatively (Signature II high) signed speckle protein gene PC1 rotations:
```python getConsistentSpeckleGenes.py speckleProteinGenes_rotations_allStudiesWithRNA_noNA.txt > speckleProteinGene_rotationChargeCounts_22cancers.txt```
Using the above file, I generated lists of consistent Signature I ("sigI_speckleProteinGenes_22cancers.txt") and Signature II high ("sigII_speckleProteinGenes_22cancers.txt") speckle protein genes that were found to have the same PC1 rotation sign in 22/22 cancer types. These speckle signature protein genes were used to generate speckle scores. 

# Generate speckle scores and create heatmaps of the speckle signature genes for each cancer type
Generates heatmaps of the speckle genes from "speckleProteinGenes_22cancers.txt", with an annotation column for the speckle score, which is calculated by the following formula:

sum((z-score sigI gene)*1/(number Sig I genes)) + sum((z-score sigII gene)*-1/(number Sig II genes)

This means that patients with Signature scores that are highly positive are highly Signature I, and patients with Signature scores that are highly negative are highly Signature II

```Rscript makeHeatmapsFromDownloadedCBioPortal.R```
The above script returns:
1. "cBioPortal_signatureHeatmaps_22cancers" - A folder containing PDFs of the speckle signature heatmaps for each speckle protein gene
2. "cBioPortal_signatureSpeckleScores_22cancers" - A folder containing the tables of speckle signature scores for each sample within each study
