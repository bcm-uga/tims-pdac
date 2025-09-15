# DNA methylation and immune infiltration mediate the impact of tobacco exposure on pancreatic adenocarcinoma outcome: a high-dimensional mediation analysis

This repository contains all scripts and resources to reproduce the analyses and figures from the manuscript:

> *DNA methylation and immune infiltration mediate the impact of tobacco exposure on pancreatic adenocarcinoma outcome: a high-dimensional mediation analysis*  

---

## Repository structure

### Root files
- **`master_simulation.R`** — runs the full simulation study.  
- **`master_pdac.R`** — runs the real data analysis (TCGA-PAAD).  
- **`master_figures.md`** — documentation and scripts to reproduce all manuscript figures.  

### Folders
- **`simulations/`** — contains scripts and functions used in the simulation framework.  
- **`real_data/`** — contains scripts for preprocessing and analyzing TCGA-PAAD data.  
- **`Figures/`** — contains scripts for generating each figure from the manuscript.  

---

## Analyses

- **Simulation study**:  
  Run with `master_simulation.R`. Evaluates the performance of HDMAX2-surv under multiple scenarios.  

- **Real data analysis (TCGA-PAAD)**:  
  Run with `master_pdac.R`. Applies HDMAX2-surv to study how tobacco exposure influences pancreatic adenocarcinoma survival via DNA methylation and immune infiltration.  

Due to size constraints, full processed results can be downloaded from **Zenodo**: [DOI: XXX].

---

## Requirements

The scripts are written in **R**.  
We recommend using R (≥ 4.2.0) with the following packages installed:  
- `survival`  
- `mediation`
- `multimediate`
- `tidyverse`  
- `hdmax2`  
- `InstaPrism`
- `MuSiC`
- `SingleCellExperiment`
- `survival`
- `ggplot2`
- `see`
- `dplyr`
- `UpSetR`
- `epimedtools`
- `reshape2`
- `mpmi`
- `reticulate`
- `ggpubr`
- `dagitty`
- `broom`
- `pheatmap`
- `survminer`
- `gridExtra`
- `foreach`
- `boot`
- `scales`
- `stringr`
- `tidyr`
- `gridExtra`
- `viridis`
- `ggrepel`
- `cowplot`
- `doBy`
- `LightLogR`
- `ComplexHeatmap`
- `patchwork`
- `glmnet`
- `HDMT`
- `MASS`
- `timereg`
- `mvtnorm`
  

---

## Usage

Clone this repository:  

```bash
git clone https://github.com/your-username/tims-pdac.git
cd tims-pdac
