################
### Figure 5 ###
################

library(ggplot2)
library(dplyr)
library(scales)

### Load and process data

# Load the serial mediation results

res_mediation1 = readRDS("real_data/results/06_TCGA_CAT1_serial_med_V2_K8_corrected_with_selected_LF.rds")
res_mediation2 = readRDS("real_data/results/06_TCGA_CAT2_serial_med_V2_K8_corrected_with_selected_LF.rds")
res_mediation = rbind(res_mediation1, res_mediation2)
dim(res_mediation)

med1 <- sapply(strsplit(res_mediation$feat, "_"), function(x) x[1])
med2 <- sapply(strsplit(res_mediation$feat, "_"), function(x) paste(x[-1], collapse = " "))

res_mediation$med1 <- med1
res_mediation$med2 <- med2

res_mediation$med2[res_mediation$med2 == "all"] <- "Tot. Imm."

idx = paste(res_mediation$med1, res_mediation$med2, sep = "-")
res_mediation$pairs = idx

DAGvalidated_CAT1 <- readRDS("real_data/results/05_DAG_validated_CAT1_withLF.rds")
DAGvalidated_CAT1 <- gsub("_", " ", DAGvalidated_CAT1)

DAGvalidated_CAT2 <- readRDS("real_data/results/05_DAG_validated_CAT2_withLF.rds")
DAGvalidated_CAT2 <- gsub("_", " ", DAGvalidated_CAT2)


### Define plot function

plot_all_pairs_all_effects <- function(res_df,
                                       sigma    = NULL,
                                       base     = 10,
                                       n_breaks = 7) {
  
  # 1) Build the data.frame without filtering
  df <- res_df %>%
    mutate(pair = paste(med1, med2, sep = " – ")) %>%
    mutate(pair = factor(pair, levels = unique(pair))) %>%
    mutate(
      sign_dir    = estimate > 0,
      significant = ifelse(ci_low > 0 | ci_high < 0, "Signif", "NS")
    )
  
  # 2) Automatically compute sigma if not provided
  if (is.null(sigma)) {
    sigma <- median(abs(df$estimate[df$estimate != 0]), na.rm = TRUE)
  }
  
  # 3) Compute pretty breaks and select 2 on each side of 0
  all_vals <- c(df$estimate, df$ci_low, df$ci_high)
  br_raw   <- pretty(all_vals, n = n_breaks)
  neg_raw  <- sort(br_raw[br_raw < 0], decreasing = TRUE)
  pos_raw  <- sort(br_raw[br_raw > 0], decreasing = FALSE)
  neg2     <- head(neg_raw, 2)
  pos2     <- head(pos_raw, 2)
  my_breaks <- c(neg2, 0, pos2)
  
  # 4) Build the plot
  p <- ggplot(df, aes(x = estimate, y = pair)) +
    # Vertical reference line at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    
    # Error bars colored by sign direction
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = sign_dir),
                   height = 0.25) +
    
    # Points colored by direction and shaped by significance
    geom_point(aes(color = sign_dir, shape = significant),
               size = 2.5) +
    
    # Pseudo-log scale with 2 negatives, 0, and 2 positives
    scale_x_continuous(
      trans  = pseudo_log_trans(sigma = sigma, base = base),
      breaks = my_breaks,
      labels = label_comma(accuracy = 1, big.mark = " "),
      guide  = guide_axis(n.dodge = 1)  # all labels on a single line
    ) +
    
    # Manual color palette
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "skyblue"),
                       guide = "none") +
    
    # Shape mapping for significance
    scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
    
    # Facet by effect type
    facet_grid(. ~ effect_type, scales = "free_y", space = "free_y", drop = FALSE) +
    
    # Titles and themes
    labs(
      title = "Effects for AMR–cell type pairs",
      x     = paste0("ACME (pseudo-log10 transform - linear scale around 0)"),
      y     = "AMR – Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y   = element_text(size = 7),
      axis.text.x   = element_text(
        size  = 7,
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.title    = element_text(size = 14),
      strip.text    = element_text(size = 11, face = "bold"),
      plot.title    = element_text(size = 15, face = "bold"),
      legend.title  = element_blank(),
      legend.text   = element_text(size = 10)
    )
  
  return(p)
}


### PANEL A: Handmade graph

### PANEL B: DAG 1

res_significatif = res_mediation[which(res_mediation$pairs %in% DAGvalidated_CAT1), ]
p_sig_pairs <- plot_all_pairs_all_effects(res_significatif)

ggsave("figures/fig5_panelB.pdf",
       p_sig_pairs, width = 11, height = 5, units = "in")
	   
### PANEL C: DAG 2


res_significatif = res_mediation[which(res_mediation$pairs %in% DAGvalidated_CAT2), ]
p_sig_pairs <- plot_all_pairs_all_effects(res_significatif)

ggsave("figures/fig5_panelC.pdf",
       p_sig_pairs, width = 11, height = 5, units = "in")
