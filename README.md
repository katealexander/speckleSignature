# Speckle Signature
This repository uses RNA-seq data to estimate nuclear speckle phenotypes. It accompanies the manuscript, Alexander et al., which is under preparation (as of 05/31/2023). Nuclear speckles are emerging as a critical layer of gene regulation. However, we lack methods to robustly define speckle phenotypes. The speckle signature is a tool that uses RNA-seq data of speckle-resident protein genes to estimate speckle phenotypes. 

# Requirements
R packages: [cBioPortalData](https://bioconductor.org/packages/release/bioc/html/cBioPortalData.html), [AnVIL](https://bioconductor.org/packages/release/bioc/html/AnVIL.html), [pheatmap](https://CRAN.R-project.org/package=pheatmap)
Python 2.7 -- The Python scripts herein are relatively simple. I expect they could easily be converted to [Python3](https://python2to3.com/), but have not tried this.

# Analysis 
### Speckle_signature_definition
Speckle_signature_definition describes how I arrived at the current working definition of the speckle signature.

### Calculate_speckle_signature_score
Calculate_speckle_signature_score describes how to generate speckle scores from RNA-seq data.

# Considerations and Limitations
The speckle signature represents one way of measuring speckle variation and does not capture all possible ways in which nuclear speckles may vary. It was generated based on consistent speckle protein gene expression patterns observed among 24 cancer types in The Cancer Genome Atlas data. However, 6 other cancer types showed distict speckle protein gene expression patterns, which are not captured by the speckle signature. Because speckle signature scores are calculated using z-scores of a group of samples, accurate calculation relies on having a large enough cohort with adequate samples displaying either speckle signature. The speckle signature score represents a relative rather than absolute qunatification of speckle signature.  

# Do RNA-seq based speckle signatures reflect speckle phenotypes?
The speckle scores calculated from RNA data matched well with imaging-based measurements of nuclear speckles, giving confidence that speckle signature scores reflect bonafide speckle phenotypes. 






