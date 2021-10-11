# PGC-PTSD-Longitudinal-Analysis
The following scripts are used to run the PGC PTSD EWAS Working Group's longitudinal analysis. The analysis steps are as follows:

## Step 1: Run Model 1 to evaluate CpGs associated with PTSS
Use the script **Model1.R** to perform a linear mixed model with DNAm values at both time points as the dependent variable, PTSS at both time points as a main effect, and a random intercept for subject, including Age,cell proportions, and ancestry PCs derived from DNAm data were included as covariates.

## Step 2: Run Model 2 to evaluate CpGs associated with *change* in PTSS
Use the script **Model2.R** to conduct a longitudinal analysis using a linear regression model, where post-deployment DNAm was modeled as a function of change in PTSS while adjusting for pre-deployment DNAm, PCs for ancestry, and changes in age (i.e. time passed between pre- and post-deployment data collection), CD8T+, CD4T+, NK, B cell, and monocyte cell proportions.

## Step 3: Meta-Analysis
Use the script **MetaAnalysis.R** to perform a meta-analysis using weighted sum of z-scores method for Model 1 and Model 2.

## Step 4: DMR Analysis
Use the script **DMRcate_usingSumStats.R** to conduct DMR analysis, using summary statisticts from each meta-analysis.


