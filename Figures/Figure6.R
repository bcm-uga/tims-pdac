################
### Figure 6 ###
################

library(ggplot2)
library(patchwork)
### Load and process data

# Get methylation data

tcga_data = readRDS("results/01_tcga_data_expo_deconv.rds")
tcga_dna = t(tcga_data$M)
smoking = ifelse(tcga_data$tobacco==0, 'Non-smoker', 'Smoker')
names(smoking) = rownames(tcga_data$M)

tcga_data$probes_info = lapply(tcga_data$probes_info, function(x) {
  names(x) = colnames(tcga_data$M)
  x
})


# Get significant AMR data

tobacco_AMR <- readRDS("results/03_tcga_significative_from_top50_tobacco_AMR_fdr0_05_V2_K8_corrected.rds")

interesting_AMR = c("AMR17", "AMR42", "AMR46")
tcga_dna_amr = lapply(interesting_AMR, function(x)
  tcga_dna[tobacco_AMR$CpG_related_info[[x]],])
names(tcga_dna_amr) = interesting_AMR

#Load AMR data
datAMR = read.csv2("results/03_tcga_AMR_mean_meth_top50_fdr0_05_V2_K8_corrected.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datAMR[] <- lapply(datAMR, function(x) as.numeric(as.character(x)))

# Load the immune cell-type estimation
datIMM = read.csv2("results/03_tcga_consensus_deconv_immune_cells.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datIMM[] <- lapply(datIMM, function(x) as.numeric(as.character(x)))
datIMM$all = rowSums(datIMM)
colnames(datIMM) = c("Macrophages", "B cells", "T cells", "NK", "DCs", "Tot. Imm." )


# Load platform
platform = readRDS("results/platform_meth_icgc.rds")
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
tcga_exprs0 = readRDS("results/study_TCGA-PAAD_trscr.rds")

genes <- c("NAA11", "CHL1", "SPTBN2")
expr = t(tcga_exprs0$data[genes, rownames(tcga_data$M)])
df = data.frame(expr,status = ifelse(tcga_data$tobacco==0, 'Non-smoker', 'Smoker'))

# Format data 

df_tcga_dna_amr = mapply(function(x,y,name) data.frame("Pos"=x,
                                                       "Group"=rep(smoking, each=nrow(y)),
                                                       "Patient"=rep(colnames(y), each=nrow(y)),
                                                       "Methylation"=c(y),
                                                       "AMR" = name ),
                         pos, tcga_dna_amr,names(pos), SIMPLIFY = F)


### Define plot functions 

plot_AMR <- function(dat, name, vec, gene) {
  
  # Plot methylation across CpG sites in an AMR region
  plot_met_AMR <- function(dat, name) {
    ggplot(dat, aes(x = Pos, y = Methylation, color = Group)) +
      geom_line(aes(group = Patient), alpha = 0.1) +
      geom_point() +
      stat_summary(
        fun = mean,
        fun.min = ~ mean(.) - sd(.),
        fun.max = ~ mean(.) + sd(.),
        geom = "ribbon",
        aes(fill = Group),
        alpha = 0.2,
        color = NA
      ) +
      stat_summary(fun = mean, geom = "line", aes(group = Group), linewidth = 1) +
      scale_color_manual(values = c("Non-smoker" = "blue", "Smoker" = "red")) +
      scale_fill_manual(values = c("Non-smoker" = "blue", "Smoker" = "red")) +
      labs(
        title = name,
        x = "CpG Position",
        y = "Methylation"
      ) +
      ylim(0, 1) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
  }
  
  # Plot correlation between AMR methylation and immune cell abundance
  plot_imm_by_AMR <- function(vec) {
    df <- do.call(rbind, lapply(vec, function(pair) {
      amr <- sub("-.*", "", pair)
      cell <- sub(".*-", "", pair)
      if (!(amr %in% colnames(datAMR)) || !(cell %in% colnames(datIMM))) return(NULL)
      
      samples <- intersect(rownames(datAMR), rownames(datIMM))
      data.frame(
        Sample = samples,
        AMR = amr,
        Cell_Type = cell,
        AMR_Value = datAMR[samples, amr],
        Cell_Prop = datIMM[samples, cell],
        Smoking_Status = smoking[samples],
        Pair = paste(amr, cell, sep = " - "),
        stringsAsFactors = FALSE
      )
    }))
    
    ggplot(df, aes(x = AMR_Value, y = Cell_Prop, color = Smoking_Status)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_smooth(method = "lm", se = FALSE) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_color_manual(values = c("Non-smoker" = "blue", "Smoker" = "red")) +
      labs(
        title = unique(df$Pair),
        x = "Mean AMR methylation",
        y = "Immune cell proportion"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  # Plot gene expression stratified by smoking status
  plot_gene_expr <- function(gene, df) {
    expr <- t(tcga_exprs0$data[genes, rownames(tcga_data$M)])
    df <- data.frame(expr, status = ifelse(tcga_data$tobacco == 0, "Non-smoker", "Smoker"))
    pval <- wilcox.test(df[[gene]] ~ df$status)$p.value
    
    ggplot(df, aes(x = status, y = df[[gene]], fill = status)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
      scale_fill_manual(values = c("Smoker" = "red", "Non-smoker" = "blue")) +
      labs(
        title = paste("Expression of", gene),
        subtitle = paste("Wilcoxon p =", signif(pval, 3)),
        x = "Smoking Status", y = "Gene Expression"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  # Generate all plots
  p_AMR <- plot_met_AMR(dat, name)
  p_imm <- plot_imm_by_AMR(vec)
  p_gene <- plot_gene_expr(gene, df)
  
  # Combine the plots into a single figure
  final_plot <- (p_AMR | p_imm | p_gene) +
    plot_layout(guides = "collect", heights = c(2, 1)) +
    plot_annotation(title = paste0("Methylation, immune infiltration, and gene expression for ", name)) &
    theme(legend.position = "bottom")
  
  return(final_plot)
}


### PANEL A: AMR17 - DCs - NAA11

name = "AMR17"
dat = df_tcga_dna_amr[["AMR17"]]
vec = c("AMR17-DCs")
gene = "NAA11"
final_plot = plot_AMR(dat, name, vec, gene)
ggsave("figures/fig6_panelA.pdf", final_plot, width = 8, height = 4, units = "in")

### PANEL B: AMR42 - Tot. Imm. - "SPTBN2"

name = "AMR42"
dat = df_tcga_dna_amr[["AMR42"]]
vec = c("AMR42-Tot. Imm.")
gene = "SPTBN2"
final_plot = plot_AMR(dat, name, vec, gene)
ggsave("figures/fig6_panelB.pdf", final_plot, width = 8, height = 4, units = "in")

### PANEL C: AMR62 - Tot. Imm. - "CHL1"

name = "AMR46"
dat = df_tcga_dna_amr[["AMR46"]]
vec = c("AMR46-Tot. Imm.")
gene = "CHL1"
final_plot = plot_AMR(dat, name, vec, gene)
ggsave("figures/fig6_panelC.pdf", final_plot, width = 8, height = 4, units = "in")
