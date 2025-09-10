library(tidyverse)
library(ggrepel)
library(cowplot)
library(doBy)
library(ggplot2)

############################################################################
## CpG mediation results
############################################################################
# load data
tcga_data = readRDS("tcga_pdac_mediation/results/01_tcga_data_expo_deconv.rds")

###################################
## acme plot  FDR 0.05
###################################

med_tobacco_dnam_CpG = readRDS("tcga_pdac_mediation/results/03_tcga_top50_tobacco_CPG_fdr0_05_V2_K8_corrected.rds")

cpg = data.frame(med_tobacco_dnam_CpG$step2_res_4)

p_ACME_005 <- ggplot(cpg, aes(
  x = est,
  y = feat
  ,
  color = est <= 0
)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5), height = 0.3) +
  geom_point(size = 2.8) +
  theme_bw() +
  labs(
    title = "ACME of FDR 5% significative CpG",
    y = "Mediators",
    x = "Est. effect [CI 2.5–97.5]"
  ) +
  theme(
    panel.border = element_blank(),
    panel.spacing = unit(0.01, "lines"),
    axis.ticks = element_blank()
  ) +
  scale_color_manual(values = c("skyblue", "red"), guide = "none")

ggsave(
  "figures/03_TCGA_tobacco_CPG_fdr0_05_acme_plot_V2_K8_corrected.pdf",
  p_ACME_005,
  width = 4,
  height = 6,
  units = "in"
)


#####################################
## acme plot FDR 0.1
####################################

med_tobacco_dnam_CpG = readRDS("results/03_tcga_top50_tobacco_CPG_fdr0_1_V2_K8_corrected.rds")

cpg = data.frame(med_tobacco_dnam_CpG$step2_res_6)

p_ACME_01 <- ggplot(cpg, aes(
  x = est,
  y = feat
  ,
  color = est <= 0
)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5), height = 0.3) +
  geom_point(size = 2.8) +
  theme_bw() +
  labs(
    title = "ACME of FDR 10% significative CpG",
    y = "Mediators",
    x = "Est. effect [CI 2.5–97.5]"
  ) +
  theme(
    panel.border = element_blank(),
    panel.spacing = unit(0.01, "lines"),
    axis.ticks = element_blank()
  ) +
  scale_color_manual(values = c("skyblue", "red"), guide = "none")

ggsave(
  "figures/03_TCGA_tobacco_CPG_fdr0_1_acme_plot_V2_K8_corrected.pdf",
  p_ACME_01,
  width = 4,
  height = 6,
  units = "in"
)
