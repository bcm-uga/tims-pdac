################
### Figure 3 ###
################


library(tidyverse)
library(ggrepel)
library(cowplot)
library(doBy)
library(ggplot2)
library(scales)  
library(LightLogR) 

# load data
med_tobacco_dnam_AMR = readRDS("real_data/results/03_tcga_top50_tobacco_AMR_fdr0_05_V2_K8_corrected.rds")
med_tobacco_dnam_CpG = readRDS("real_data/results/02_tcga_med_tobacco_dnam_V2_K8_corrected.rds")
pf = readRDS("real_data/results/03_platform_info_tcga_meth.rds")

tcga_data = readRDS("real_data/results/01_tcga_data_expo_deconv.rds")

############################################################################
## manhattan plot
############################################################################
# prepare data

max2_pvalues = med_tobacco_dnam_CpG$hdmax2_step1_param$max2_pvalues
cpg = data.frame(rownames(pf), pf$Chromosome, pf$Start, max2_pvalues)
colnames(cpg) = c("cpg", "chrom", "pos", "pval")
cpg$chrom <- as.integer(gsub("chr", "", cpg$chrom))
cpg <- orderBy( ~ chrom, cpg)
cpg$index <- 1:nrow(cpg)



pf$Gene_Symbol_unique = unname(sapply(pf$Gene_Symbol, function(x) {
  paste(unique(unlist(strsplit(x, split = ";"))), collapse = ";")
}))
annot_df <- pf
annot_df$gene <- pf$Gene_Symbol_unique


get_major_genes_for_amrs <- function(amr_list, annot_df) {
  # Initialize the result vector
  major_genes <- character(length(amr_list))
  
  # Loop over each AMR
  for (i in seq_along(amr_list)) {
    probes <- amr_list[[i]]
    
    # Select genes associated with probes of this AMR
    matched_genes <- annot_df$gene[rownames(annot_df)%in% probes]
    
    # Clean: split multiple names if separated by ";"
    genes_split <- unlist(strsplit(matched_genes, split = ";"))
    genes_clean <- genes_split[genes_split != ""]
    
    if (length(genes_clean) == 0) {
      major_genes[i] <- NA
    } else {
      gene_tab <- table(genes_clean)
      major_genes[i] <- names(gene_tab)[which.max(gene_tab)]
    }
  }
  
  return(major_genes)
}

amr_list = med_tobacco_dnam_AMR$cpg_related_50
# Result: a vector with 1 gene per AMR
gene_amr <- get_major_genes_for_amrs(amr_list, annot_df)


cpg_related = amr_list

# retrieve first CpG of the amr to get a position
for (i in 1:length(cpg_related)) {
  cpg_related[[i]] = cpg_related[[i]][1]
}

amr = data_frame(est = med_tobacco_dnam_AMR$step2_res_50$est,
                 cpg = cpg_related,
                 pval_amr = med_tobacco_dnam_AMR$AMR_info$p)

amr$index = cpg$index[which(cpg$cpg %in% amr$cpg)]
amr$mark = "AMR"
amr$gene = gene_amr

# Clean and convert chromosome format
cpg$chrom <- gsub("chr", "", cpg$chrom)
cpg$chrom <- as.numeric(cpg$chrom)
# Sort
cpg <- cpg[order(cpg$chrom, cpg$pos), ]

# Get unique sorted chromosomes
chroms <- sort(unique(cpg$chrom))

# Maximum size per chromosome
chr_max <- tapply(cpg$pos, cpg$chrom, max)

# Calculate offsets
chr_offset <- numeric(length(chroms))
names(chr_offset) <- chroms

for (i in 2:length(chroms)) {
  chr_offset[i] <- chr_offset[i - 1] + chr_max[as.character(chroms[i - 1])]
}


cpg$index <- cpg$pos + chr_offset[as.character(cpg$chrom)]
amr$index <- cpg$index[match(amr$cpg, cpg$cpg)]
chr_labels <- tapply(cpg$index, cpg$chrom, function(x) median(range(x)))

pos_chr <- NULL
for (i in 1:length(table(cpg$chrom))) {
  cpg2 <- cpg[cpg$chrom == i, ]
  pos_chr <- rbind(pos_chr, cbind(chrom = i, index = median(cpg2$index)))
}
pos_chr <- as.data.frame(pos_chr)

chr_1 <- cpg[cpg$chrom %in% seq(1, 21, 2), ]
chr_2 <- cpg[cpg$chrom %in% seq(2, 22, 2), ]

# level harmonisation
cpg$mark <- factor("CpG", levels = c("CpG", "AMR" ))
chr_1$mark <- factor("CpG", levels = c("CpG", "AMR"))
chr_2$mark <- factor("CpG", levels = c("CpG", "AMR"))
amr$mark <- factor("AMR", levels = c("CpG", "AMR"))
pos_chr$mark <- factor("CpG", levels = c("CpG", "AMR"))
amr$amr_id = med_tobacco_dnam_AMR$step2_res_50$feat

## plot

p_man = ggplot() +
  geom_point(data = chr_1, aes(index, -log10(pval)), col = "grey80") +
  geom_point(data = chr_2, aes(index, -log10(pval)), col = "grey70") +
  geom_text(
    data = pos_chr,
    aes(x = index, y = 0.5, label = chrom),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  geom_point(data = amr, aes(index, -log10(pval_amr)), color = "slategray4", size = 3)+
  geom_segment(data = amr, aes(x = index, xend = index, y = 0, yend = -log10(pval_amr)), color = "slategray4", size = 0.6)+
  geom_text_repel(data = amr, aes(index, -log10(pval_amr), label = amr_id),
                  color = "black", fontface = "italic") +
  facet_wrap( ~ mark, scales = "free_y", nrow = 2) +
  xlab("Chromosome position") +
  ylab("-log P") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
  )
p_man
# save plot
ggsave("figures/03_TCGA_tobacco_AMR_fdr0_05_man_V2_K8_corrected.pdf", p_man, width = 12, height = 5, units = "in")

####################################
## acme plot
####################################

amr$acme_pval <- med_tobacco_dnam_AMR$step2_res_50$pval
amr$CI_2.5 <- med_tobacco_dnam_AMR$step2_res_50$CI_2.5
amr$CI_97.5 <- med_tobacco_dnam_AMR$step2_res_50$CI_97.5
amr$signif <- ifelse(amr$acme_pval <= 0.05, "pval <= 0.05", "n.s.")
amr$amr_id = med_tobacco_dnam_AMR$step2_res_50$feat
amr$surv_dist = med_tobacco_dnam_AMR$all_step2_res$effects$univariate$best_distribution

ordered_amr_feat <- amr$amr_id[order(amr$est)]

# pseudo log version

p_ACME_log <- ggplot(amr, aes(
  x     = est,
  y     = factor(amr_id, levels = ordered_amr_feat),
  color = surv_dist,
  #color = est <= 0,
  shape = signif
)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5), height = 0.3) +
  geom_point(size = 2.8) +
  theme_bw() +
  labs(
    title = "A. ACME (Average Causal Mediation Effect)",
    y     = "Mediators",
    x     = "Effect (pseudo-log10)"
  ) +
  scale_color_manual(
    values = c("magenta3", "mediumseagreen" , "orange2"),
    name   = "Survival distribution"
  )+
  scale_shape_manual(
    values = c("pval <= 0.05" = 16, "n.s." = 17),
    guide  = guide_legend(title = "Significance")
  ) +
  scale_x_continuous(
    trans   = pseudo_log_trans(base = 10),
    breaks  = c(-100, -10, -1, 0, 1, 10, 100),
    labels  = function(x) x
  )

ggsave(
  "figures/03_TCGA_tobacco_AMR_acme_pseudolog.pdf",
  p_ACME_log,
  width  = 6,
  height = 8,
  units  = "in"
)



####################################################
# Indirect effect decomposition
######################################################


amr$XM_est = med_tobacco_dnam_AMR$all_step2_res$effects$univariate$xm$Estimate
amr$MY_est = med_tobacco_dnam_AMR$all_step2_res$effects$univariate$my$Estimate


box <- data.frame(name = c("Smoking", "Survival"),
                  X = c(3.6,6.4),  # before it was 0 and 10
                  Y = c(0, 0))

med <- subset(amr, acme_pval <= 0.05)
med <- med[order(-med$est), ] 

med$X <- 5
med$Y <- seq(10, -10, length.out = nrow(med))

# from tobacco to genes
med$X_XM_seg_s <- 3.6
med$Y_XM_seg_s <- 0
med$X_XM_seg_e <- 4.8
med$Y_XM_seg_e <- med$Y

# from genes to survival
med$X_MY_seg_s <- 5.2
med$Y_MY_seg_s <- med$Y
med$X_MY_seg_e <- 6.4
med$Y_MY_seg_e <- 0

# Colors and styles
med$fill_acme <- ifelse(med$est <= 0, "Negative", "Positive")
med$col_xm <- ifelse(med$XM_est <= 0, "Negative", "Positive")
med$col_my <- ifelse(med$MY_est <= 0, "Negative", "Positive")

# plot
plot_decompo = ggplot() +
  geom_segment(
    data = med,
    aes(
      x = X_XM_seg_s,
      y = Y_XM_seg_s,
      xend = X_XM_seg_e,
      yend = Y_XM_seg_e,
      col = col_xm
    ),
    size = 1
  ) +
  geom_segment(
    data = med,
    aes(
      x = X_MY_seg_s,
      y = Y_MY_seg_s,
      xend = X_MY_seg_e,
      yend = Y_MY_seg_e,
      col = col_my
    ),
    size = 1
  ) +
  geom_label(
    data = box,
    aes(X, Y, label = name),
    size = 7,
    fontface = "bold"
  ) +
  geom_label(
    data = med,
    aes(X, Y, label = amr_id, fill = fill_acme),
    fontface = "italic",
    size = 5
  ) +
  scale_color_manual(values = c("red", "skyblue")) +
  scale_fill_manual(values = c("mistyrose", "aliceblue")) +
  labs(color = "Regression coeff sign", fill = "Indirect Effect (ACME)") +
  coord_cartesian(clip = "off") +
  expand_limits(x = c(3.3, 6.7)) + 
  theme_void() +
  theme(
    legend.position = c(0.05, 0.95),
    # top left
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 10)
  )

ggsave("figures/03_ACME_effect_decomposition.pdf",plot_decompo, height = 10, width = 7)
