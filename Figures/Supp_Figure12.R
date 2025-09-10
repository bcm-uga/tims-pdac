######################
### Supp Figure 12 ###
######################

library(ggplot2)
library(gridExtra)

### Load and process data

# Get methylation data

tcga_data = readRDS("tcga_pdac_mediation/results/01_tcga_data_expo_deconv.rds")
tcga_dna = t(tcga_data$M)
smoking = ifelse(tcga_data$tobacco==0, 'Non-smoker', 'Smoker')
names(smoking) = rownames(tcga_data$M)

tcga_data$probes_info = lapply(tcga_data$probes_info, function(x) {
  names(x) = colnames(tcga_data$M)
  x
})


# Get significant AMR data

tobacco_AMR <- readRDS("tcga_pdac_mediation/results/03_tcga_significative_from_top50_tobacco_AMR_fdr0_05_V2_K8_corrected.rds")

interesting_AMR = c("AMR1", "AMR17", "AMR42", "AMR46")
tcga_dna_amr = lapply(interesting_AMR, function(x)
  tcga_dna[tobacco_AMR$CpG_related_info[[x]],])
names(tcga_dna_amr) = interesting_AMR

#Load AMR data
datAMR = read.csv2("tcga_pdac_mediation/results/03_tcga_AMR_mean_meth_top50_fdr0_05_V2_K8_corrected.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datAMR[] <- lapply(datAMR, function(x) as.numeric(as.character(x)))


# Load platform
platform = readRDS("tcga_pdac_mediation/data/TCGA_PAAD/platform_meth_icgc.rds")
cpg_annotation = platform

cpg = lapply(tcga_dna_amr, rownames)
length(unlist(cpg))
cpg_info_table = lapply(cpg, function(x) {
  res = cpg_annotation[x,]
  res = res[!duplicated(res$Start),]
  res = res[!duplicated(res$End),]
  res
})

pos = lapply(cpg, function(x) cpg_annotation[x,]$Start)
chr = sapply(cpg, function(x) unique(cpg_annotation[x,]$Chromosome))
feat_type = mapply(function(x, y) cpg_annotation[cpg_annotation$Start>=min(x) & cpg_annotation$End<=(max(x)+1) & cpg_annotation$Chromosome==y, c("Chromosome","Start","End","Feature_Type")],
                   pos, chr, SIMPLIFY = F)

#load gene expression
tcga_exprs0 = readRDS("tcga_pdac_mediation/data/TCGA_PAAD/study_TCGA-PAAD_trscr.rds")

genes <- c("HOXC4", "TSC2", "FERMT3", "PIK3R1") # Known genes of interest
expr = t(tcga_exprs0$data[genes, rownames(tcga_data$M)])
df = data.frame(expr,status = ifelse(tcga_data$tobacco==0, 'Non-smoker', 'Smoker'))


### Plot Figure


cols <- c("Smoker" = "red", "Non-smoker" = "blue")

plot_list <- list()

for (gene in genes) {
  # Wilcoxon test
  pval <- wilcox.test(df[[gene]] ~ df$status)$p.value
  
  # Plot
  p <- ggplot(df, aes_string(x = "status", y = gene, fill = "status")) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
    scale_fill_manual(values = cols) +
    labs(
      title = paste("Expression of", gene),
      subtitle = paste("Wilcoxon p =", signif(pval, 3)),
      x = "Smoking Status", y = "Gene Expression"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
  
  plot_list[[gene]] <- p
}

pdf("figures/supp_fig12.pdf", width = 10, height = 12)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()
