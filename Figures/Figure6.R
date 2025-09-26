################
### Figure 6 ###
################

library(ggplot2)
library(patchwork)
### Load and process data

# Get methylation data

tcga_data = readRDS("real_data/results/01_tcga_data_expo_deconv.rds")
tcga_dna = t(tcga_data$M)
smoking = tcga_data$tobacco
names(smoking) = rownames(tcga_data$M)

tcga_data$probes_info = lapply(tcga_data$probes_info, function(x) {
  names(x) = colnames(tcga_data$M)
  x
})


# Get significant AMR data

tobacco_AMR <- readRDS("real_data/results/03_tcga_significative_from_top50_tobacco_AMR_fdr0_05_V2_K8_corrected.rds")

interesting_AMR = c("AMR4", "AMR5", "AMR9", "AMR17", "AMR20","AMR23", "AMR27")
tcga_dna_amr = lapply(interesting_AMR, function(x)
  tcga_dna[tobacco_AMR$CpG_related_info[[x]],])
names(tcga_dna_amr) = interesting_AMR

#Load AMR data
datAMR = read.csv2("real_data/results/03_tcga_AMR_mean_meth_top50_fdr0_05_V2_K8_corrected.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datAMR[] <- lapply(datAMR, function(x) as.numeric(as.character(x)))

# Load the immune cell-type estimation
datIMM = read.csv2("real_data/results/03_tcga_consensus_deconv_immune_cells.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datIMM[] <- lapply(datIMM, function(x) as.numeric(as.character(x)))
datIMM$all = rowSums(datIMM)
colnames(datIMM) = c("Macrophages", "B cells", "T cells", "NK", "DCs", "Tot. Imm." )


# Load platform
platform = readRDS("real_data/tcga_pdac_mediation/data/platform_meth_icgc.rds")
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
tcga_exprs0 = readRDS("real_data/tcga_pdac_mediation/data/TCGA_PAAD/study_TCGA-PAAD_trscr.rds")

genes <- c("ZBTB42", "KCNQ1", "CLGN", "PIWIL1", "RAD51AP2", "FAT4", "PRDM16") #genes linked to AMRs
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
						       scale_color_manual(
						         values = c(
						           "Lifelong Non-Smoker" = "#08519c",                  
						           "Current Reformed Smoker for > 15 yrs" = "#6baed6",
						           "Current Reformed Smoker for < or = 15 yrs" = "#fb6a4a",
						           "Current Smoker" = "#a50f15"
						         ),
						         labels = c(
						           "Lifelong Non-Smoker" = "Non-Smoker",
						           "Current Reformed Smoker for > 15 yrs" = "Reformed >15 yrs",
						           "Current Reformed Smoker for < or = 15 yrs" = "Reformed ≤15 yrs",
						           "Current Smoker" = "Smoker"
						         ),
						         breaks = c(
						           "Lifelong Non-Smoker",
						           "Current Reformed Smoker for > 15 yrs",
						           "Current Reformed Smoker for < or = 15 yrs",
						           "Current Smoker"
						         )
						       ) +
						       scale_fill_manual(
						         values = c(
						           "Lifelong Non-Smoker" = "#08519c",                 
						           "Current Reformed Smoker for > 15 yrs" = "#6baed6",
						           "Current Reformed Smoker for < or = 15 yrs" = "#fb6a4a",
						           "Current Smoker" = "#a50f15"
						         ),
						         labels = c(
						           "Lifelong Non-Smoker" = "Non-Smoker",
						           "Current Reformed Smoker for > 15 yrs" = "Reformed >15 yrs",
						           "Current Reformed Smoker for < or = 15 yrs" = "Reformed ≤15 yrs",
						           "Current Smoker" = "Smoker"
						         ),
						         breaks = c(
						           "Lifelong Non-Smoker",
						           "Current Reformed Smoker for > 15 yrs",
						           "Current Reformed Smoker for < or = 15 yrs",
						           "Current Smoker"
						         )
						       ) +
						       labs(
						         title = name,
						         x = "CpG Position",
						         y = "Methylation"
						       ) +
						       ylim(0, 1) +
						       theme_minimal() +
						       theme(
						         axis.text.x = element_text(angle = 45, hjust = 1),
						         legend.position = "bottom",
						         legend.title = element_blank(),
						         legend.box = "horizontal",
						         legend.justification = "center"
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
						       geom_smooth(method = "lm", se = FALSE, color = "black") + 
						       scale_x_continuous(limits = c(0, 1)) +
						       scale_color_manual(
						         values = c(
						           "Lifelong Non-Smoker" = "#08519c",                  
						           "Current Reformed Smoker for > 15 yrs" = "#6baed6",
						           "Current Reformed Smoker for < or = 15 yrs" = "#fb6a4a",
						           "Current Smoker" = "#a50f15"
						         ),
						         labels = c(
						           "Lifelong Non-Smoker" = "Non-Smoker",
						           "Current Reformed Smoker for > 15 yrs" = "Reformed >15 yrs",
						           "Current Reformed Smoker for < or = 15 yrs" = "Reformed ≤15 yrs",
						           "Current Smoker" = "Smoker"
						         ),
						         breaks = c(
						           "Lifelong Non-Smoker",
						           "Current Reformed Smoker for > 15 yrs",
						           "Current Reformed Smoker for < or = 15 yrs",
						           "Current Smoker"
						         )
						       ) +
						       labs(
						         title = unique(df$Pair),
						         x = "Mean AMR methylation",
						         y = "Immune cell proportion"
						       ) +
						       theme_minimal() +
						       theme(
						         legend.position = "none",
						         legend.title = element_blank(),
						         legend.box = "horizontal",
						         legend.justification = "center"
						       )
						   }
  
  
						   # Plot gene expression stratified by smoking status
						   plot_gene_expr <- function(gene, df) {
						     library(ggpubr)
    
						     
						     expr <- t(tcga_exprs0$data[genes, rownames(tcga_data$M)])
						     df <- data.frame(expr, status = tcga_data$tobacco)
    
						    
						     df$status <- factor(df$status,
						                         levels = c(
						                           "Lifelong Non-Smoker",
						                           "Current Reformed Smoker for > 15 yrs",
						                           "Current Reformed Smoker for < or = 15 yrs",
						                           "Current Smoker"
						                         ),
						                         labels = c("Non-Smoker", "Reformed >15 yrs", "Reformed ≤15 yrs", "Smoker")
						     )
    
						    
						     comparisons <- list(
						       c("Non-Smoker", "Reformed >15 yrs"),
						       c("Non-Smoker", "Reformed ≤15 yrs"),
						       c("Non-Smoker", "Smoker"),
						       c("Reformed >15 yrs", "Reformed ≤15 yrs"),
						       c("Reformed >15 yrs", "Smoker"),
						       c("Reformed ≤15 yrs", "Smoker")
						     )
    
						     ggplot(df, aes(x = status, y = .data[[gene]], fill = status)) +
						       geom_violin(trim = FALSE, alpha = 0.5) +
						       geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
						       scale_fill_manual(
						         values = c(
						           "Non-Smoker" = "#08519c",
						           "Reformed >15 yrs" = "#6baed6",
						           "Reformed ≤15 yrs" = "#fb6a4a",
						           "Smoker" = "#a50f15"
						         )
						       ) +
						       labs(
						         title = paste("Expression of", gene),
						         x = "Smoking Status",
						         y = "Gene Expression"
						       ) +
						       theme_minimal() +
						       theme(
						         legend.position = "none",  # <-- légende supprimée
						         axis.text.x = element_text(angle = 45, hjust = 1)
						       ) +
						       stat_compare_means(comparisons = comparisons, method = "wilcox.test")
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

						 ### PANEL A: 

						 name = "AMR4"
						 dat = df_tcga_dna_amr[["AMR4"]]
						 vec = c("AMR4-B cells")
						 gene = "ZBTB42"
						 final_plot = plot_AMR(dat, name, vec, gene)
						 ggsave("figures/fig6_panelA.pdf", final_plot, width = 8, height = 6, units = "in")

						 ### PANEL B: 

						 name = "AMR9"
						 dat = df_tcga_dna_amr[["AMR9"]]
						 vec = c("AMR9-Tot. Imm.")
						 gene = "CLGN"
						 final_plot = plot_AMR(dat, name, vec, gene)
						 ggsave("figures/fig6_panelB.pdf", final_plot, width = 8, height = 6, units = "in")

						 ### PANEL C:

						 name = "AMR20"
						 dat = df_tcga_dna_amr[["AMR20"]]
						 vec = c("AMR20-Tot. Imm.")
						 gene = "RAD51AP2"
						 final_plot = plot_AMR(dat, name, vec, gene)
						 ggsave("figures/fig6_panelC.pdf", final_plot, width = 8, height = 6, units = "in")

						 ### PANEL D: 

						 name = "AMR27"
						 dat = df_tcga_dna_amr[["AMR27"]]
						 vec = c("AMR27-Tot. Imm.")
						 gene = "PRDM16"
						 final_plot = plot_AMR(dat, name, vec, gene)
						 ggsave("figures/fig6_panelD.pdf", final_plot, width = 8, height = 6, units = "in")
