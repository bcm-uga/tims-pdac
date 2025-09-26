# Required packages
library(survival)
library(boot)
library(dplyr)

# # Main function: Serial mediation with survival outcome
# serial_mediation_survie2 <- function(expo, mediateur1, mediateur2, covars, time, status, R = 1000) {
#   
#   # Prepare analysis dataset
#   data_to_analyse <- data.frame(
#     expo = expo,
#     m1 = mediateur1,
#     m2 = mediateur2,
#     time = as.numeric(time),
#     status = as.numeric(status),
#     covars
#   )
#   
#   # Remove individuals with time = 0
#   data_to_analyse <- data_to_analyse[data_to_analyse$time > 0, ]
#   data_to_analyse$id <- 1:nrow(data_to_analyse)
#   
#   # Function to estimate effects
#   theta <- function(data, index) {
#     dat <- data[index, ]
#     
#     # Models for mediators
#     fitM1 <- lm(m1 ~ expo + ., data = dat[, c("m1", "expo", colnames(covars))])
#     fitM2 <- lm(m2 ~ expo + m1 + ., data = dat[, c("m2", "expo", "m1", colnames(covars))])
#     
#     # Automatically select best survival distribution using AIC
#     dists <- c("weibull", "exponential", "loglogistic", "lognormal")
#     aics <- sapply(dists, function(d) {
#       tryCatch({
#         AIC(survreg(Surv(time, status) ~ expo + m1 + m2 + ., 
#                     data = dat[, c("time", "status", "expo", "m1", "m2", colnames(covars))],
#                     dist = d))
#       }, error = function(e) Inf)
#     })
#     best_dist <- dists[which.min(aics)]
#     
#     # Survival model with selected distribution
#     fitY <- survreg(Surv(time, status) ~ expo + m1 + m2 + ., 
#                     data = dat[, c("time", "status", "expo", "m1", "m2", colnames(covars))], 
#                     dist = best_dist)
#     
#     # Create counterfactual dataset (4 scenarios of a1, a2, a3)
#     dat_ext <- dat[rep(1:nrow(dat), each = 4), ]
#     dat_ext$a1 <- rep(rep(0:1, each = 2), times = nrow(dat))
#     dat_ext$a2 <- rep(0:1, times = 2 * nrow(dat))
#     dat_ext$a3 <- rep(0:1, each = 2 * nrow(dat))
#     
#     # Predict mediators
#     dat_ext$pred_m1_a1 <- predict(fitM1, newdata = within(dat_ext, { expo = a1 }))
#     dat_ext$pred_m2_a2m1a1 <- predict(fitM2, newdata = within(dat_ext, {
#       expo = a2
#       m1 = pred_m1_a1
#     }))
#     
#     # Predict survival outcome
#     dat_ext$Y_hat <- predict(fitY, newdata = within(dat_ext, {
#       expo = a3
#       m1 = pred_m1_a1
#       m2 = pred_m2_a2m1a1
#     }), type = "response")
#     
#     # Estimate effects
#     fitNEM <- lm(Y_hat ~ a3 + a1 + a2, data = dat_ext)
#     coef_fit <- coef(fitNEM)
#     
#     total_effect <- coef_fit["a3"] + coef_fit["a1"] + coef_fit["a2"]
#     direct_effect <- coef_fit["a3"]
#     indirect_m1 <- coef_fit["a1"]
#     indirect_m2 <- coef_fit["a2"]
#     indirect_joint <- indirect_m1 + indirect_m2
#     
#     return(c(direct_effect, indirect_m1, indirect_m2, indirect_joint, total_effect))
#   }
#   
#   # Bootstrap estimation
#   set.seed(42)
#   boot_res <- boot(data = data_to_analyse, statistic = theta, R = R)
#   
#   linfunCI.perc <- function(boot.out, index, conf = 0.90) {
#     est <- boot.out$t0[index]
#     alpha <- (1 - conf) / 2
#     v <- boot.out$t[, index]
#     ci <- quantile(v, probs = c(alpha, 1 - alpha), na.rm = TRUE)
#     c(estimate = est, ci_low = ci[1], ci_high = ci[2])
#   }
#   
#   
#   names_effets <- c("Direct effect", "Indirect effect via M1", "Indirect effect via M2",
#                     "Total indirect effect", "Total effect")
#   results <- t(sapply(1:5, function(i) linfunCI.perc(boot_res, i)))
#   rownames(results) <- names_effets
#   
#   return(round(results, 6))
# }
# 
# # Function to apply mediation to all DMR × Cell type pairs
# run_serial_mediation_grid2 <- function(
#     mat_dmr,
#     mat_cell,
#     expo,
#     covars,
#     time,
#     status,
#     R = 1000
# ) {
#   results_list <- list()
#   
#   for (i in 1:ncol(mat_dmr)) {
#     for (j in 1:ncol(mat_cell)) {
#       
#       cat("Running mediation:", colnames(mat_dmr)[i], "x", colnames(mat_cell)[j], "\n")
#       
#       res <- tryCatch({
#         serial_mediation_survie2(
#           expo = expo,
#           mediateur1 = mat_dmr[, i],
#           mediateur2 = mat_cell[, j],
#           covars = covars,
#           time = time,
#           status = status,
#           R = R
#         )
#       }, error = function(e) {
#         message("Error for", colnames(mat_dmr)[i], "x", colnames(mat_cell)[j], ":", e$message)
#         return(NULL)
#       })
#       
#       if (!is.null(res)) {
#         res_df <- as.data.frame(res)
#         res_df$med1 <- colnames(mat_dmr)[i]
#         res_df$med2 <- colnames(mat_cell)[j]
#         res_df$effect_type <- rownames(res)
#         results_list[[length(results_list) + 1]] <- res_df
#       }
#     }
#   }
#   
#   if (length(results_list) > 0) {
#     all_results <- do.call(rbind, results_list)
#     rownames(all_results) <- NULL
#     return(all_results)
#   } else {
#     warning("No valid result generated.")
#     return(NULL)
#   }
# }
# 


# Function to apply mediation to all DMR × Cell type pairs

run_serial_mediation_grid2_alter <- function(
    mat_dmr,
    mat_cell,
    expo,
    covars_LF,
    covars,
    time,
    status,
    R = 1000
) {
  results_list <- list()
  
  for (i in 1:ncol(mat_dmr)) {
    for (j in 1:ncol(mat_cell)) {
      
      cat("Running mediation:", colnames(mat_dmr)[i], "x", colnames(mat_cell)[j], "\n")
      
      res <- tryCatch({
        serial_mediation_survie2_alter(
          expo = expo,
          mediateur1 = mat_dmr[, i],
          mediateur2 = mat_cell[, j],
          covars_LF = covars_LF,
          covars = covars,
          time = time,
          status = status,
          R = R
        )
      }, error = function(e) {
        message("Error for", colnames(mat_dmr)[i], "x", colnames(mat_cell)[j], ":", e$message)
        return(NULL)
      })
      
      if (!is.null(res)) {
        res_df <- as.data.frame(res)
        res_df$med1 <- colnames(mat_dmr)[i]
        res_df$med2 <- colnames(mat_cell)[j]
        res_df$effect_type <- rownames(res)
        results_list[[length(results_list) + 1]] <- res_df
      }
    }
  }
  
  if (length(results_list) > 0) {
    all_results <- do.call(rbind, results_list)
    rownames(all_results) <- NULL
    return(all_results)
  } else {
    warning("No valid result generated.")
    return(NULL)
  }
}


plot_effects_for_celltype <- function(res_df, cell_type, effet = "Effet indirect via M2") {
  
  # Filtrer les résultats pour un type cellulaire donné et un effet donné
  df_plot <- subset(res_df, med2 == cell_type & effect_type == effet)
  
  # Extraire le numéro de DMR (par ex. "DMR12" → 12)
  df_plot$dmr_num <- as.numeric(gsub("AMR", "", df_plot$med1))
  
  # Trier par numéro de DMR
  df_plot <- df_plot[order(df_plot$dmr_num), ]
  
  # Pour que les AMR s'affichent dans le bon ordre sur l'axe y
  df_plot$med1 <- factor(df_plot$med1, levels = df_plot$med1)
  
  # Créer le plot
  p <- ggplot(df_plot, aes(x = estimate, y = med1)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3) +
    geom_point(size = 2.5, aes(color = estimate > 0)) +
    scale_color_manual(values = c("skyblue", "red"), guide = "none") +
    labs(
      title = paste("Mediated effect via", cell_type),
      x = "Estimated effect[IC 95%]",
      y = "AMR"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 15, face = "bold")
    )
  
  return(p)
}

plot_effects_for_celltype_tot <- function(res_df, cell_type, effet = "Effet indirect total") {
  
  # Filtrer les résultats pour un type cellulaire donné et un effet donné
  df_plot <- subset(res_df, med2 == cell_type & effect_type == effet)
  
  # Extraire le numéro de DMR (par ex. "DMR12" → 12)
  df_plot$dmr_num <- as.numeric(gsub("DMR", "", df_plot$med1))
  
  # Trier par numéro de DMR
  df_plot <- df_plot[order(df_plot$dmr_num), ]
  
  # Pour que les AMR s'affichent dans le bon ordre sur l'axe y
  df_plot$med1 <- factor(df_plot$med1, levels = df_plot$med1)
  
  # Créer le plot
  p <- ggplot(df_plot, aes(x = estimate, y = med1)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3) +
    geom_point(size = 2.5, aes(color = estimate > 0)) +
    scale_color_manual(values = c("skyblue", "red"), guide = "none") +
    labs(
      title = paste("Effets médiés via", cell_type),
      x = "Effet estimé [IC 95%]",
      y = "DMR"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 15, face = "bold")
    )
  
  return(p)
}



# plot_facet_all_effects <- function(res_df, effet_types = NULL, cell_types = NULL) {
#   
#   # Filtrer si nécessaire
#   if (!is.null(effet_types)) {
#     res_df <- res_df %>% filter(effect_type %in% effet_types)
#   }
#   if (!is.null(cell_types)) {
#     res_df <- res_df %>% filter(med2 %in% cell_types)
#   }
#   
#   # Extraire le numéro de DMR
#   res_df$dmr_num <- as.numeric(gsub("AMR", "", res_df$med1))
#   
#   # Trier par DMR
#   res_df <- res_df %>% arrange(dmr_num)
#   
#   # Facteur pour l'axe Y (ordre des DMR)
#   res_df$med1 <- factor(res_df$med1, levels = unique(res_df$med1))
#   
#   # Créer variable de significativité : IC n'inclut pas 0
#   res_df$significant <- ifelse(res_df$ci_low > 0 | res_df$ci_high < 0, "Signif", "NS")
#   
#   # Créer le plot
#   p <- ggplot(res_df, aes(x = estimate, y = med1)) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
#     geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
#     geom_point(aes(color = estimate > 0, shape = significant), size = 2.5) +
#     scale_color_manual(values = c("skyblue", "red"), guide = "none") +
#     scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
#     facet_grid(effect_type ~ med2, scales = "free_y", space = "free_y") +
#     labs(
#       title = "Mediated effects according to cell type and AMR",
#       x = "Estimated effect [CI 95%]",
#       y = "AMR"
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.y = element_text(size = 8),
#       axis.text.x = element_text(size = 10),
#       axis.title = element_text(size = 14),
#       strip.text = element_text(size = 11, face = "bold"),
#       plot.title = element_text(size = 15, face = "bold"),
#       legend.title = element_text(size = 12),
#       legend.text = element_text(size = 10)
#     )
#   
#   return(p)
# }

library(ggplot2)
library(dplyr)
library(scales)   # pour pseudo_log_trans()

plot_facet_all_effects <- function(res_df,
                                   effet_types = NULL,
                                   cell_types  = NULL,
                                   sigma       = NULL,
                                   base        = 10,
                                   breaks      = NULL) {
  
  # 0) Calcul automatique de sigma si non fourni
  #    (par exemple la médiane des |estimate| non-nuls)
  if (is.null(sigma)) {
    sigma <- median(abs(res_df$estimate[res_df$estimate != 0]), na.rm = TRUE)
  }
  
  # 1) Filtrage
  if (!is.null(effet_types)) {
    res_df <- res_df %>% filter(effect_type %in% effet_types)
  }
  if (!is.null(cell_types)) {
    res_df <- res_df %>% filter(med2 %in% cell_types)
  }
  
  # 2) Extraction du numéro de DMR et factorage
  res_df <- res_df %>%
    mutate(
      dmr_num     = as.numeric(gsub("AMR", "", med1)),
      med1        = factor(med1, levels = unique(med1)),
      significant = ifelse(ci_low > 0 | ci_high < 0, "Signif", "NS")
    ) %>%
    arrange(dmr_num)
  
  # 3) Calcul des breaks si non fournis : on prend des jolis « pretty »
  if (is.null(breaks)) {
    all_vals <- c(res_df$estimate, res_df$ci_low, res_df$ci_high)
    mag      <- pretty(all_vals, n = 7)
    # on s’assure d’avoir zéro et ±
    breaks   <- unique(c(mag[mag < 0], 0, mag[mag > 0]))
  }
  
  # 4) Construction du plot
  p <- ggplot(res_df, aes(x = estimate, y = med1)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
    geom_point(aes(color = estimate > 0, shape = significant),
               size = 2.5) +
    
    # Pseudo-log sur l’axe x
    scale_x_continuous(
      trans  = pseudo_log_trans(base = base, sigma = sigma),
      breaks = breaks,
      labels = function(x) x
    ) +
    
    scale_color_manual(values = c("skyblue", "red"), guide = "none") +
    scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
    
    facet_grid(effect_type ~ med2, scales = "free_y", space = "free_y") +
    
    labs(
      title = "Mediated effects selon type de cellule et AMR",
      x     = sprintf("Estimated effect (pseudo-log10, σ = %.2g)", sigma),
      y     = "AMR"
    ) +
    theme_bw() +
    theme(
      axis.text.y   = element_text(size = 8),
      axis.text.x   = element_text(size = 10),
      axis.title    = element_text(size = 14),
      strip.text    = element_text(size = 11, face = "bold"),
      plot.title    = element_text(size = 15, face = "bold"),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 10)
    )
  
  return(p)
}




library(ggplot2)
library(dplyr)
library(scales)

plot_all_effects_as_facets <- function(res_df,
                                       sigma     = NULL,
                                       base      = 10,
                                       n_breaks  = 7) {
  # 1) Ne garder que med2 == "all"
  df <- res_df %>%
    filter(med2 == "all") %>%
    mutate(
      # Extraction du numéro pour trier les AMR
      dmr_num     = as.numeric(gsub("AMR", "", med1)),
      med1        = factor(med1, levels = unique(med1)),
      significant = ifelse(ci_low > 0 | ci_high < 0, "Signif", "NS")
    ) %>%
    arrange(dmr_num)
  
  # 2) Calcul automatique de sigma si besoin (médiane des |estimate|)
  if (is.null(sigma)) {
    sigma <- median(abs(df$estimate[df$estimate != 0]), na.rm = TRUE)
  }
  
  # 3) Détermination des breaks “jolis” sur l’échelle d’origine
  all_vals   <- c(df$estimate, df$ci_low, df$ci_high)
  pretty_raw <- pretty(all_vals, n = n_breaks)
  breaks     <- unique(c(pretty_raw[pretty_raw < 0], 0, pretty_raw[pretty_raw > 0]))
  
  # 4) Construction du plot
  ggplot(df, aes(x = estimate, y = med1)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
    geom_point(aes(color = estimate > 0, shape = significant), size = 2.5) +
    
    # Pseudo-log (compressions extrêmes / zone linéaire) mais labels en unités d'origine
    scale_x_continuous(
      trans  = scales::pseudo_log_trans(sigma = sigma, base = base),
      breaks = breaks,
      labels = function(x) x
    ) +
    scale_color_manual(values = c("skyblue", "red"), guide = "none") +
    scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
    
    # On remplace cell_types par effect_type en colonnes
    facet_grid(. ~ effect_type, scales = "free_y", space = "free_y") +
    
    labs(
      title = "Tous les effets pour med2 = all",
      x     = sprintf("Estimate (pseudo-log10, σ = %.2g)", sigma),
      y     = "AMR"
    ) +
    theme_bw() +
    theme(
      axis.text.y   = element_text(size = 8),
      axis.text.x   = element_text(size = 10),
      axis.title    = element_text(size = 14),
      strip.text    = element_text(size = 11, face = "bold"),
      plot.title    = element_text(size = 15, face = "bold"),
      legend.title  = element_text(size = 12),
      legend.text   = element_text(size = 10)
    )
}


library(dplyr)
library(ggplot2)
library(scales)

plot_signif_pairs_all_effects <- function(res_df,
                                          sigma    = NULL,
                                          base     = 10,
                                          n_breaks = 7) {
  # 1) Détecter les paires (med1, med2) où
  #    effect_type == "Total indirect effect" & IC95 ne croise pas 0
  sel_pairs <- res_df %>%
    filter(
      effect_type == "Total indirect effect",
      ci_low     > 0 | ci_high < 0
    ) %>%
    transmute(pair = paste(med1, med2, sep = " – ")) %>%
    distinct() %>%
    pull(pair)
  
  # 2) Construire la table filtrée et ordonner le factor 'pair'
  df <- res_df %>%
    mutate(pair = paste(med1, med2, sep = " – ")) %>%
    filter(pair %in% sel_pairs) %>%
    # On ordonne 'pair' dans l'ordre d'apparition
    mutate(pair = factor(pair, levels = unique(pair))) %>%
    # Significativité sur chaque effet
    mutate(significant = ifelse(ci_low > 0 | ci_high < 0, "Signif", "NS"))
  
  # 3) Calcul automatique de sigma si non fourni
  if (is.null(sigma)) {
    sigma <- median(abs(df$estimate[df$estimate != 0]), na.rm = TRUE)
  }
  
  # 4) Déterminer des bornes “pretty” pour l'axe X en valeurs d'origine
  all_vals   <- c(df$estimate, df$ci_low, df$ci_high)
  br_raw     <- pretty(all_vals, n = n_breaks)
  breaks     <- unique(c(br_raw[br_raw < 0], 0, br_raw[br_raw > 0]))
  
  # 5) Construction du plot
  p <- ggplot(df, aes(x = estimate, y = pair)) +
    # Ligne verticale à 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    # Barres d'erreur horizontales
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
    # Points colorés suivant le signe et la significativité
    geom_point(aes(color = estimate > 0, shape = significant),
               size = 2.5) +
    
    # Échelle pseudo-log
    scale_x_continuous(
      trans  = scales::pseudo_log_trans(sigma = sigma, base = base),
      # 1) Choix de breaks « pretty » réduit
      breaks = scales::breaks_pretty(n = 5),
      # 2) label_comma pour séparer les milliers
      labels = scales::label_comma(accuracy = 1, big.mark = " "),
      # 3) guide_axis pour répartir sur 2 rangées
      guide  = guide_axis(n.dodge = 2)
    ) +
    scale_color_manual(values = c("skyblue", "red"), guide = "none") +
    scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
    
    # Facettes : une colonne par type d'effet
    facet_grid(. ~ effect_type, scales = "free_y", space = "free_y", drop = FALSE) +
    
    # Titres et thèmes
    labs(
      title = "Effects for significant AMR–cell type pairs",
      x     = sprintf("Estimate (pseudo-log10, σ = %.2g)", sigma),
      y     = "AMR – Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y   = element_text(size = 8),
      axis.text.x   = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
      axis.title    = element_text(size = 14),
      strip.text    = element_text(size = 11, face = "bold"),
      plot.title    = element_text(size = 15, face = "bold"),
      legend.title  = element_blank(),
      legend.text   = element_text(size = 10)
    )
  
  return(p)
}


plot_all_pairs_all_effects <- function(res_df,
                                       sigma    = NULL,
                                       base     = 10,
                                       n_breaks = 7) {
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  # 1) Construction du data.frame sans filtrage
  df <- res_df %>%
    mutate(pair = paste(med1, med2, sep = " – ")) %>%
    mutate(pair = factor(pair, levels = unique(pair))) %>%
    mutate(
      sign_dir    = estimate > 0,
      significant = ifelse(ci_low > 0 | ci_high < 0, "Signif", "NS")
    )
  
  # 2) Calcul automatique de sigma si non fourni
  if (is.null(sigma)) {
    sigma <- median(abs(df$estimate[df$estimate != 0]), na.rm = TRUE)
  }
  
  # 3) Détermination des "pretty" breaks et sélection 2 de part et d’autre + 0
  all_vals <- c(df$estimate, df$ci_low, df$ci_high)
  br_raw   <- pretty(all_vals, n = n_breaks)
  neg_raw  <- sort(br_raw[br_raw < 0], decreasing = TRUE)
  pos_raw  <- sort(br_raw[br_raw > 0], decreasing = FALSE)
  neg2     <- head(neg_raw, 2)
  pos2     <- head(pos_raw, 2)
  my_breaks <- c(neg2, 0, pos2)
  
  # 4) Construction du plot
  p <- ggplot(df, aes(x = estimate, y = pair)) +
    # Ligne verticale à 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    
    # Barres d'erreur colorées selon sign_dir
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = sign_dir),
                   height = 0.25) +
    
    # Points colorés selon sign_dir et formes selon significant
    geom_point(aes(color = sign_dir, shape = significant),
               size = 2.5) +
    
    # Échelle pseudo-log avec 2 négatifs, 0, 2 positifs
    scale_x_continuous(
      trans  = pseudo_log_trans(sigma = sigma, base = base),
      breaks = my_breaks,
      labels = label_comma(accuracy = 1, big.mark = " "),
      guide  = guide_axis(n.dodge = 1)  # même ligne pour tous les labels
    ) +
    
    # Palette d'origine
    scale_color_manual(values = c("FALSE" = "red", "TRUE" =  "skyblue"),
                       guide = "none") +
    
    # Formes pour significativité
    scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
    
    # Facettes par type d'effet
    facet_grid(. ~ effect_type, scales = "free_y", space = "free_y", drop = FALSE) +
    
    # Titres et thèmes
    labs(
      title = "Effects for AMR–cell type pairs",
      x     = sprintf("Estimate (pseudo-log10, σ = %.2g)", sigma),
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

##################################
# Required packages
library(survival)
library(boot)
library(dplyr)

# Main function: Serial mediation with survival outcome
serial_mediation_survie2_alter <- function(expo, mediateur1, mediateur2, covars_LF, covars, time, status, R = 1000) {
  
  # Prepare analysis dataset
  data_to_analyse <- data.frame(
    expo = expo,
    m1 = mediateur1,
    m2 = mediateur2,
    time = as.numeric(time),
    status = as.numeric(status),
    covars_LF,
    covars
  )
  
  # Remove individuals with time = 0
  data_to_analyse <- data_to_analyse[data_to_analyse$time > 0, ]
  data_to_analyse$id <- 1:nrow(data_to_analyse)
  
  # Function to estimate effects
  theta <- function(data, index) {
    dat <- data[index, ]
    
    # Models for mediators
    fitM1 <- lm(m1 ~ expo + ., data = dat[, c("m1", "expo", colnames(covars))])
    fitM2 <- lm(m2 ~ expo + m1 + ., data = dat[, c("m2", "expo", "m1", colnames(covars_LF))])
    
    # Automatically select best survival distribution using AIC
    dists <- c("weibull", "exponential", "loglogistic", "lognormal")
    aics <- sapply(dists, function(d) {
      tryCatch({
        AIC(survreg(Surv(time, status) ~ expo + m1 + m2 + ., 
                    data = dat[, c("time", "status", "expo", "m1", "m2", colnames(covars))],
                    dist = d))
      }, error = function(e) Inf)
    })
    best_dist <- dists[which.min(aics)]
    
    # Survival model with selected distribution
    fitY <- survreg(Surv(time, status) ~ expo + m1 + m2 + ., 
                    data = dat[, c("time", "status", "expo", "m1", "m2", colnames(covars))], 
                    dist = best_dist)
    
    # Create counterfactual dataset (4 scenarios of a1, a2, a3)
    dat_ext <- dat[rep(1:nrow(dat), each = 4), ]
    dat_ext$a1 <- rep(rep(0:1, each = 2), times = nrow(dat))
    dat_ext$a2 <- rep(0:1, times = 2 * nrow(dat))
    dat_ext$a3 <- rep(0:1, each = 2 * nrow(dat))
    
    # Predict mediators
    dat_ext$pred_m1_a1 <- predict(fitM1, newdata = within(dat_ext, { expo = a1 }))
    dat_ext$pred_m2_a2m1a1 <- predict(fitM2, newdata = within(dat_ext, {
      expo = a2
      m1 = pred_m1_a1
    }))
    
    # Predict survival outcome
    dat_ext$Y_hat <- predict(fitY, newdata = within(dat_ext, {
      expo = a3
      m1 = pred_m1_a1
      m2 = pred_m2_a2m1a1
    }), type = "response")
    
    # Estimate effects
    fitNEM <- lm(Y_hat ~ a3 + a1 + a2, data = dat_ext)
    coef_fit <- coef(fitNEM)
    
    total_effect <- coef_fit["a3"] + coef_fit["a1"] + coef_fit["a2"]
    direct_effect <- coef_fit["a3"]
    indirect_m1 <- coef_fit["a1"]
    indirect_m2 <- coef_fit["a2"]
    indirect_joint <- indirect_m1 + indirect_m2
    
    return(c(direct_effect, indirect_m1, indirect_m2, indirect_joint, total_effect))
  }
  
  # Bootstrap estimation
  set.seed(42)
  boot_res <- boot(data = data_to_analyse, statistic = theta, R = R)
  
  # Compute bootstrap percentile confidence intervals
  linfunCI.perc <- function(boot.out, index, conf = 0.90) {
    est <- boot.out$t0[index]
    alpha <- (1 - conf) / 2
    v <- boot.out$t[, index]
    ci <- quantile(v, probs = c(alpha, 1 - alpha), na.rm = TRUE)
    c(estimate = est, ci_low = ci[1], ci_high = ci[2])
  }
  
  names_effets <- c("Direct effect", "Indirect effect via M1", "Indirect effect via M2",
                    "Total indirect effect", "Total effect")
  results <- t(sapply(1:5, function(i) linfunCI.perc(boot_res, i)))
  rownames(results) <- names_effets
  
  return(round(results, 6))
}

# Function to apply mediation to all DMR × Cell type pairs

run_serial_mediation_grid2_alter <- function(
    mat_dmr,
    mat_cell,
    expo,
    covars_LF,
    covars,
    time,
    status,
    R = 1000
) {
  results_list <- list()
  
  for (i in 1:ncol(mat_dmr)) {
    for (j in 1:ncol(mat_cell)) {
      
      cat("Running mediation:", colnames(mat_dmr)[i], "x", colnames(mat_cell)[j], "\n")
      
      res <- tryCatch({
        serial_mediation_survie2_alter(
          expo = expo,
          mediateur1 = mat_dmr[, i],
          mediateur2 = mat_cell[, j],
          covars_LF = covars_LF,
          covars = covars,
          time = time,
          status = status,
          R = R
        )
      }, error = function(e) {
        message("Error for", colnames(mat_dmr)[i], "x", colnames(mat_cell)[j], ":", e$message)
        return(NULL)
      })
      
      if (!is.null(res)) {
        res_df <- as.data.frame(res)
        res_df$med1 <- colnames(mat_dmr)[i]
        res_df$med2 <- colnames(mat_cell)[j]
        res_df$effect_type <- rownames(res)
        results_list[[length(results_list) + 1]] <- res_df
      }
    }
  }
  
  if (length(results_list) > 0) {
    all_results <- do.call(rbind, results_list)
    rownames(all_results) <- NULL
    return(all_results)
  } else {
    warning("No valid result generated.")
    return(NULL)
  }
}


plot_all_pairs_all_effects <- function(res_df,
                                       sigma    = NULL,
                                       base     = 10,
                                       n_breaks = 7) {
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(stringr)
  
  # 1) Préparation du data.frame (feat contient déjà les noms)
  df <- res_df %>%
    mutate(
      amr_num = as.numeric(stringr::str_extract(feat, "(?<=AMR)\\d+")),
      sign_dir = estimate > 0,
      significant = ifelse(ci_low > 0 | ci_high < 0, "Signif", "NS")
    ) %>%
    arrange(amr_num) %>%
    mutate(feat = factor(feat, levels = unique(feat)))
  
  # 2) Calcul automatique de sigma si non fourni
  if (is.null(sigma)) {
    sigma <- median(abs(df$estimate[df$estimate != 0]), na.rm = TRUE)
  }
  
  # 3) Détermination des "pretty" breaks
  all_vals <- c(df$estimate, df$ci_low, df$ci_high)
  br_raw   <- pretty(all_vals, n = n_breaks)
  neg_raw  <- sort(br_raw[br_raw < 0], decreasing = TRUE)
  pos_raw  <- sort(br_raw[br_raw > 0], decreasing = FALSE)
  neg2     <- head(neg_raw, 2)
  pos2     <- head(pos_raw, 2)
  my_breaks <- c(neg2, 0, pos2)
  
  # 4) Construction du plot
  p <- ggplot(df, aes(x = estimate, y = feat)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = sign_dir),
                   height = 0.25) +
    geom_point(aes(color = sign_dir, shape = significant),
               size = 2.5) +
    scale_x_continuous(
      trans  = pseudo_log_trans(sigma = sigma, base = base),
      breaks = my_breaks,
      labels = label_comma(accuracy = 1, big.mark = " "),
      guide  = guide_axis(n.dodge = 1)
    ) +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" =  "skyblue"),
                       guide = "none") +
    scale_shape_manual(values = c("Signif" = 16, "NS" = 17)) +
    facet_grid(. ~ effect_type, scales = "free_y", space = "free_y", drop = FALSE) +
    labs(
      title = "Effects for AMR–cell type pairs",
      x     = sprintf("Estimate (pseudo-log10, σ = %.2g)", sigma),
      y     = "AMR – Cell type"
    ) +
    theme_bw() +
    theme(
      axis.text.y   = element_text(size = 7),
      axis.text.x   = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
      axis.title    = element_text(size = 14),
      strip.text    = element_text(size = 11, face = "bold"),
      plot.title    = element_text(size = 15, face = "bold"),
      legend.title  = element_blank(),
      legend.text   = element_text(size = 10)
    )
  
  return(p)
}
