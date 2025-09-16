#################
### tims-pdac ###
################
# Master script to run real data application of TGCA-PAAD
# to run this pipeline create a results/ folder in existing real_data/ folder
# for 00_sc_processing.Rmd you need to load from zenodo tcga_pdac_mediation/data/genref_hadaca3/00_peng_k_2019.rds, tcga_pdac_mediation/data/genref_hadaca3/00_raghavan_s_2021.rds
# for 01_TCGA_processing_deconv.Rmd you need to load from zenodo tcga_pdac_mediation/data/cometh_lot1/transcriptome/T_raw.rds, tcga_pdac_mediation/data/TCGA_PAAD/study_TCGA-PAAD_meth.rds, tcga_pdac_mediation/data/TCGA_PAAD/study_TCGA-PAAD_trscr.rds, tcga_pdac_mediation/data/TCGA_PAAD/PublishedSampleTable.csv, tcga_pdac_mediation/data/TCGA_PAAD/clinical.project-tcga-paad.2024-07-30/exposure.tsv


# ## Preprocessing
# # Single cell data preprocessing
rmarkdown::render("real_data/00_sc_processing.Rmd", run_pandoc = FALSE, clean = TRUE)

## Step 1 Deconvolution
# We developed a consensus deconvolution method that integrates both bulk and single-cell reference datasets to enhance accuracy and robustness. 
rmarkdown::render("real_data/01_TCGA_processing_deconv.Rmd", run_pandoc = FALSE, clean = TRUE)

## Step 2 Mediation
# High-dimensional mediation analysis: identification of potential DNA methylation probes mediating the pathway between tobacco consumption and survival in PDAC patients.
rmarkdown::render("real_data/02_TCGA_mediation.Rmd", run_pandoc = FALSE, clean = TRUE)

## Step 3 AMR research and AMR mediation effects estimation
# Detection of Aggregated Methylation Regions among DNA probes p-value  and estimation of mediated AMR effects
rmarkdown::render("real_data/03_TCGA_survival_deconv_tobacco_AMR.Rmd", run_pandoc = FALSE, clean = TRUE)

## Step 4 Correlation
# Correlation study between latent factors, immune profil, AMR
rmarkdown::render("real_data/04_TCGA_immune_AMR_correlation.Rmd", run_pandoc = FALSE, clean = TRUE)

## Step 5 Causality
#We implemented a causal discovery framework based on conditional independence tests to identify causal relationships among the following variables: tobacco exposure, mediating AMR, total immune infiltration, B cell infiltration proportion, dendritic cell infiltration proportion, and survival. 
rmarkdown::render("real_data/05_causality.Rmd", run_pandoc = FALSE, clean = TRUE)

## Step 6 Serial mediation
# We sought to determine whether indirect effects exist between tobacco exposure and survival in PDAC patients through mediator chains (AMR–immune infiltration pairs identified by causal discovery), acknowledging that in the context of multiple mediators the total effect can be decomposed by considering them simultaneously.
rmarkdown::render("real_data/06_TCGA_serial_mediation_V2.Rmd", run_pandoc = FALSE, clean = TRUE)

# ## Step 7 Interpretation
#rmarkdown::render("real_data/07_AMR_biology.Rmd", run_pandoc = FALSE)
