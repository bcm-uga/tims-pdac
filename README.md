# DNA methylation and immune infiltration mediate the impact of tobacco exposure on pancreatic adenocarcinoma outcome: a high-dimensional mediation analysis

This repository contains all scripts and resources to reproduce the analyses and figures from the manuscript:

> *DNA methylation and immune infiltration mediate the impact of tobacco exposure on pancreatic adenocarcinoma outcome: a high-dimensional mediation analysis*  

---

## Repository structure

### Analysis

- **`master_simulation.R`** — runs the simulation study.  
- **`master_pdac.R`** — runs the analysis on real data (TCGA-PAAD).
- **`master_figures.md`** — details which results and scripts files are required to generate the figures.  

Due to size constraints, full processed results can be downloaded from **Zenodo**: [DOI: XXX].

---

## Requirements

The scripts are written in **R**.  
We recommend using R (≥ 4.2.0) with the following packages installed:  
- `survival`  
- `glmnet`  
- `mediation`  
- `data.table`  
- `tidyverse`  
- `igraph`  
- (and any additional packages listed at the top of each script).  

---

## Usage

Clone this repository:  

```bash
git clone https://github.com/your-username/tims-pdac.git
cd tims-pdac
