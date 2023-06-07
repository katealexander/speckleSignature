# Relating speckle signature to pateint data
This pipeline describes how I used expression of speckle signature speckle protein genes to assess cancer patient survival. 

# Download patient data
The following script downloads patient data and places into the folder "patientData_cBioPortal":

```Rscript downloadPatientData_cBioPortal.R```

This script generates the folder, "patientData_cBioPortal", which is included in this repository, downloaded 06/05/2023.

# Calculate survival based on speckle signature
The following Rscript makes use of speckle scores stored in "Speckle_signature_definition/cBioPortal_signatureSpeckleScores_22cancers" and uses the patiend data downloaded above. 

```Rscript survival_cBioPortal_fromDownloadedData.R```

The above script returns a folder "cBioPortal_survival" that contains three folders of Kaplan Meier plots:
1. "survival_IvsII" assesses survival by positive and negative speckle scores.
2. "survival_topQuarterScores" assesses survival by the top and bottom 25% of speckle scores.
3. "categoricalData_speckleNtile_bargraphs" assesses whether other categorical variables in patient data correlate with speckle signature by spliting patients into quartiles by speckle score and performing Fishers exact tests. 

Desease-specific (DSS) and overall (OS) survival are both assessed, with file names starting with "DSS" or "OS", respectively. 

# Calculate speckle signature impact on survival in early- and late-stage KIRC TCGA cohort
The following script generates Kaplan Meier survival plots for the KIRC TCGA cohort with one plot for stage I and II and one plot for stage III and IV, with patients split by speckle Signature

```Rscript kirc_survival_cBioPortal_grade.R```

These plots are stored in the folder, "cBioPortal_survival/survival_byGrade"

# Calculate speckle signature impact on survival in VHL wild type and mutant KIRC TCGA cohort
The survival analysis showed that speckle signature predicts survival in the KIRC renal carcinoma cohort. Since clear cell renal cell carcinoma, which comprises the KIRC TCGA cohort, frequently harbors inactivating copy number loss or missense mutations in the VHL protein, I assessed survival in the VHL wild type and VHL mutant KIRC cohort. 
1. Downloaded Copy Number Variant number data for VHL in the KIRC TCGA cohort from the [GDC](https://portal.gdc.cancer.gov/). Not included in repository because it's not in a very useful format
2. Used R script "make_VHL_cnv_caseFile.R" to generate "KIRC_VHL_CNV.txt", provided in this repository. Copy number loss 
3. Downloaded VHL simple nucleotide variation sample list in the KIRC TCGA cohort from the GDC, called "KIRC_TCGA_VHL_casesEffected.2023-02-24.tsv"
4. Downloaded the case sets that were assayed for simple nucleotide variation ("case_set_TCGA_KIRC__Simple_Nucleotide_Variation.2023-02-24.tsv") and copy number variation ("case_TCGA_KIRC__Copy_Number_Variation__Gene_Level_Copy_Number.tsv")
5. In the below script, samples were considered VHL mutant if they contained a simple nucleotide variation or copy number loss in VHL and wild type if they did not. 

```Rscript kirc_VHLmutant_survival_cBioportalData.R```

