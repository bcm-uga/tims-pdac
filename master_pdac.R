#################
### tims-pdac ###
#################

# Master script to run real data application of TGCA-PAAD

## Preprocessing
# Single cell data preprocessing
knitr::knit_child("real_data/00_sc_processing.Rmd")

## Step 1 Deconvolution
# We developed a consensus deconvolution method that integrates both bulk and single-cell reference datasets to enhance accuracy and robustness. 
knitr::knit_child("real_data/01_TCGA_processing_deconv.Rmd")

## Step 2 Mediation
# High-dimensional mediation analysis: identification of potential DNA methylation probes mediating the pathway between tobacco consumption and survival in PDAC patients
knitr::knit_child("real_data/02_TCGA_mediation.Rmd")

## Step 3 AMR research and AMR mediation effects estimation
# Detection of Aggregated Methylation Regions among DNA probes p-value  and estimation of mediated AMR effects
knitr::knit_child("real_data/03_TCGA_survival_deconv_tobacco_AMR.Rmd")

## Step 4 Correlation
# Correlation study between latent factors, immune profil, AMR
knitr::knit_child("real_data/04_TCGA_immune_AMR_correlation.Rmd")

## Step 5 Causality
#We implemented a causal discovery framework based on conditional independence tests to identify causal relationships among the following variables: tobacco exposure, mediating AMR, total immune infiltration, B cell infiltration proportion, dendritic cell infiltration proportion, and survival 
knitr::knit_child("real_data/05_causality.Rmd")

## Step 6 Serial mediation
# We sought to determine whether indirect effects exist between tobacco exposure and survival in PDAC patients through mediator chains (AMR–immune infiltration pairs identified by causal discovery), acknowledging that in the context of multiple mediators the total effect can be decomposed by considering them simultaneously.
knitr::knit_child("real_data/06_TCGA_serial_mediation.Rmd")

## Step 7 Interpretation
knitr::knit_child("real_data/07_AMR_biology.Rmd")
