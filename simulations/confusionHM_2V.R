
library(ggplot2)
library(stringr)
library(patchwork)
library(grid)
library(gridExtra)
#library(VennDiagram)

rdsFileNames = read.csv(input_rds_files, header = FALSE, col.names = c("step2", "mediators"))
sim = as.numeric(str_match(rdsFileNames$step2, ".*_sim_([0-9]+).*")[,2])
rdsFileNames = cbind(rdsFileNames, sim)
rdsFileNames <- rdsFileNames[order(rdsFileNames$step2),]




#########
#########
# definir param1 param2 dans rule et wildcard
#########


#construction de la liste de résultats correspondant au parametre du csv chargé (simulation, resultatmediation et numero du réplicat)
results = vector(mode = "list")
 for (f in seq_along(rdsFileNames$step2)) {
    res_effect = readRDS(rdsFileNames$step2[f])
    #res_effects[[f]] = res_effect
    sim_effect = readRDS(rdsFileNames$mediators[f])
    #sim_effects[[f]] = sim_effect
    results[[f]] = list(res_effect = res_effect ,sim_effect=sim_effect, replicat = sim[f])
    #results$sim_effects[[f]] = sim_effect
 }




##########
# A COMMENTER
############
# param1 = "mean_A"
# param2 = "mean_B"
# param3 = "n"
################



param1Name = param1
param2Name = param2



param1vars = param2vars = models = simus = c()


 for (r in seq_along(results)){
   param1var = results[[r]]$sim_effect$param_values[[param1Name]]
   param2var = results[[r]]$sim_effect$param_values[[param2Name]]
   
   param1vars = c(param1vars, param1var)
   param2vars = c(param2vars, param2var)
   
   models = c(models, results[[r]]$sim_effect$model)
   simus = c(simus, results[[r]]$replicat)
 }
   
param1vars = unique(param1vars)
param1vars = sort(param1vars)
param2vars = unique(param2vars)
param2vars = sort(param2vars)
models = unique(models)
simus = unique(simus)

# 
# select_results = function(results, param1 , param2, model, simu){
#   for(i in 1:length(results)) {
#     # if(results[[i]]$sim_effect$param_values[[param1Name]]==param1 && results[[i]]$sim_effect$param_values[[param2Name]]==param2 && results[[i]]$sim_effect$model==model && results[[i]]$replicat==simu) 
#     #   {
#       mediators = results[[i]]$sim_effect$mediators
#       res_feat = results[[i]]$res_effect$effects$univariate$ACME$feat
#       inter = intersect(names(mediators), res_feat)
#       effet_A = as.numeric(results[[i]]$sim_effect$A_effect[results[[i]]$sim_effect$mediators])
#       effet_B = as.numeric(results[[i]]$sim_effect$B_effect[results[[i]]$sim_effect$mediators])
#       sign_simu = sign(effet_A*effet_B)
#       sign_simu[sign_simu == -1] = "-"
#       sign_simu[sign_simu == 1] = "+"
#       res_med = as.numeric(results[[i]]$res_effect$effects$univariate$ACME$est)
#       sign_res = sign(res_med)
#       sign_res[sign_res == -1] = "-"
#       sign_res[sign_res == 1] = "+"
#       names_all_probes = union(names(mediators), res_feat) 
#       
#       data = data.frame(names_all_probes= names_all_probes, simulation= rep("NS", length(names_all_probes)), resultats = rep("NC", length(names_all_probes)))
#       
#       data$simulation[which(names_all_probes%in%names(mediators))]= sign_simu
#       data$resultats[which(names_all_probes%in%res_feat)]= sign_res
#       
#       return(data)
#    # }
#   }
# }


select_results <- function(results, param1Value, param2Value, model, simu) {
  for (i in seq_along(results)) {
    r <- results[[i]]

    p1_match <- r$sim_effect$param_values[[param1Name]] == param1Value
    p2_match <- r$sim_effect$param_values[[param2Name]] == param2Value
    model_match <- r$sim_effect$model == model
    simu_match <- r$replicat == simu

    if (p1_match && p2_match && model_match && simu_match) {
      mediators_name <- names(r$sim_effect$mediators)
      
      res_feat <- r$res_effect$effects$univariate$ACME$feat
      inter = intersect(mediators_name, res_feat)
      effet_A <- as.numeric(r$sim_effect$A_effect[r$sim_effect$mediators])
      effet_B <- as.numeric(r$sim_effect$B_effect[r$sim_effect$mediators])
      sign_simu <- sign(effet_A * effet_B)
      sign_simu[sign_simu == -1] <- "-"
      sign_simu[sign_simu == 1] <- "+"

      res_med <- as.numeric(r$res_effect$effects$univariate$ACME$est)
      sign_res <- sign(res_med)
      sign_res[sign_res == -1] <- "-"
      sign_res[sign_res == 1] <- "+"

      all_probes <- union(mediators_name, res_feat)

      data <- data.frame(
        names_all_probes = all_probes,
        simulation = rep("NC", length(all_probes)),
        resultats = rep("NS", length(all_probes))
      )

      data$simulation[all_probes %in% mediators_name] <- sign_simu
      data$resultats[all_probes %in% res_feat] <- sign_res

      return(data)
    }
    
  }

  # Aucun match trouvé : lancer un warning
  stop(sprintf("Aucun match trouvé pour param1 = %s, param2 = %s, model = %s, simu = %s", param1Value, param2Value, model, simu))

}

# Initialisation de la liste des graphiques
# plots <- list()
# model_colors <- c("param" = "lightblue", "aalen" = "pink")
#
#
#
# # Boucles pour créer les graphiques
#  # Alterner d'abord par modèle
#   for (param1var in param1vars) {
#     for (param2var in param2vars) {
#       for (model in models) {
#       # Créer une data frame simulée (remplacer par tes données réelles)
#       data_sign = data.frame()
#       colnames <- c("simu", "result")
#
#       for (s in 1:length(simus)) {
#         data <- select_results(results = results, param1Value = param1var, param2Value = param2var , model = model, simu = s-1)
#         data_sign = rbind(data_sign, data)
#       }
#
#       # Créer une matrice de confusion
#       confusion_matrix <- table(data_sign$simulation, data_sign$resultats)
#       plot_data <- as.data.frame(confusion_matrix)
#       colnames(plot_data) <- c("simulation", "resultats", "Freq")
#
#       # Générer le graphique
#       plot_data$Freq <- plot_data$Freq / length(simus)
#       plot_data$Freq[plot_data$simulation == "NS" & plot_data$resultats == "NC"] <- NA
#
#
#       p1 <- ggplot(plot_data, aes(x = simulation, y = resultats, fill = Freq)) +
#         geom_tile(color = "white", na.rm = TRUE) +
#         theme_bw() +
#         coord_equal() +
#         scale_fill_gradient(high = scales::muted(model_colors[model]), low = model_colors[model], na.value = "white") +
#         guides(fill = FALSE) +
#         geom_text(aes(label = round(Freq, 0)), color = "black", size = 4, na.rm = TRUE) +
#         theme(axis.text.x = element_text(size = 4),
#               axis.text.y = element_text(size = 4),
#               axis.title.x = element_text(size = 4),
#               axis.title.y = element_text(size = 4),
#               plot.title = element_text(size = 8)) +
#         ggtitle(paste("Model:", model))
#
#
#
#       # Ajouter le graphique à la liste
#       plots <- append(plots, list(p1))
#
#     }
#   }
# }
#
# # Réorganisation pour alterner par méthode
# # Récupérer les indices des modèles
# param_plots <- plots[seq(1, length(plots), by = 2)] # Tous les plots pour "param"
# aalen_plots <- plots[seq(2, length(plots), by = 2)] # Tous les plots pour "aalen"
#
# # Alterner les plots
# alternating_plots <- c(rbind(param_plots, aalen_plots))
#
# # Créer la grille de graphiques
# plot_grid <- arrangeGrob(
#   grobs = alternating_plots,
#   ncol = 6 # Nombre de colonnes pour afficher les paires
# )
#
# # Ajouter des annotations globales pour les axes
# x_label <- textGrob(paste0(param1, " ",param1vars[1],  " ",param1vars[2]), gp = gpar(fontsize = 10))
# y_label <- textGrob(paste0(param2, " ",param2vars[1],  " ",param2vars[2]), gp = gpar(fontsize = 10), rot = 90)
# # Combiner la grille des graphiques avec les annotations globales
# final_plot <- grid.arrange(
#   arrangeGrob(
#     y_label,
#     plot_grid,
#     ncol = 2,
#     widths = c(1, 40)  # Ajuster les proportions pour laisser de l'espace à l'axe Y
#   ),
#   x_label,
#   nrow = 2,
#   heights = c(40, 1)  # Ajuster les proportions pour laisser de l'espace à l'axe X
# )
#
# # Afficher le graphique final
# grid.newpage()
# grid.draw(final_plot)
#
#
#
# ##################################
# # ggsave( "plot_last.pdf" ,final_plot)
#
# ggsave(output_filename1, final_plot)
#
#

library(ggplot2)
library(gridExtra)
library(grid)

# Initialisation
plots <- list()
model_colors <- c("param" = "lightblue", "aalen" = "pink")

# Stockage dans une liste nommée : [[param2]][[param1]][[model]]
for (param2var in param2vars) {
  for (param1var in param1vars) {
    for (model in models) {
      
      data_sign <- data.frame()
      
      for (s in seq_along(simus)) {
        data <- select_results(
          results = results,
          param1Value = param1var,
          param2Value = param2var,
          model = model,
          simu = simus[s]
        )
        data_sign <- rbind(data_sign, data)
      }

      confusion_matrix <- table(data_sign$simulation, data_sign$resultats)
      plot_data <- as.data.frame(confusion_matrix)
      colnames(plot_data) <- c("simulation", "resultats", "Freq")
      plot_data$Freq <- plot_data$Freq / length(simus)
      plot_data$Freq[plot_data$simulation == "NS" & plot_data$resultats == "NC"] <- NA

      p <- ggplot(plot_data, aes(x = simulation, y = resultats, fill = Freq)) +
        geom_tile(color = "white", na.rm = TRUE) +
        scale_fill_gradient(
          high = scales::muted(model_colors[model]),
          low = model_colors[model],
          na.value = "white"
        ) +
        coord_equal() +
        geom_text(aes(label = round(Freq, 0)), size = 4, na.rm = TRUE) +
        theme_minimal(base_size = 8) +
        labs(
          x = "Simulation",
          y = "Résultats",
          title = paste("Model:", model),
          subtitle = paste(param1, "=", param1var, "|", param2, "=", param2var)
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 10, face = "bold")
        )
      
      plots[[paste(param2var, param1var, model, sep = "_")]] <- p
    }
  }
}

# Réorganiser les plots pour les empiler par lignes/colonnes
# => chaque ligne = un param2, chaque colonne = un param1
grid_plots <- list()
for (param2val in param2vars) {
  row_plots <- list()
  for (param1val in param1vars) {
    row_plots <- c(row_plots, list(
      plots[[paste(param2val, param1val, "param", sep = "_")]],
      plots[[paste(param2val, param1val, "aalen", sep = "_")]]
    ))
  }
  grid_plots <- c(grid_plots, row_plots)
}

# Créer la grille 2x2 ou plus selon le nombre de paramètres
final_plot <- arrangeGrob(
  grobs = grid_plots,
  nrow = length(param2vars),
  ncol = 2 * length(param1vars)  # 2 colonnes par valeur de param1 (param + aalen)
)

# Affichage et export
grid.newpage()
grid.draw(final_plot)
ggsave(output_filename1, final_plot, width = 10, height = 8)

#
# library(ggplot2)
# library(gridExtra)
# library(grid)
# library(stringr)
# library(dplyr)
#
# # ===================================
# # Définition des paramètres analysés
# # ===================================
#
# param1Name <- param1
# param2Name <- param2
#
# # ===================================
# # Lecture des fichiers
# # ===================================
# rdsFileNames <- read.csv(input_rds_files, header = FALSE, col.names = c("step2", "mediators"))
# rdsFileNames$sim <- as.numeric(str_match(rdsFileNames$step2, ".*_sim_([0-9]+).*")[,2])
# rdsFileNames <- rdsFileNames[order(rdsFileNames$step2),]
#
# # ===================================
# # Chargement des résultats
# # ===================================
# results <- vector("list", length(rdsFileNames$step2))
# for (f in seq_along(rdsFileNames$step2)) {
#   res_effect <- readRDS(rdsFileNames$step2[f])
#   sim_effect <- readRDS(rdsFileNames$mediators[f])
#   results[[f]] <- list(res_effect = res_effect, sim_effect = sim_effect, replicat = rdsFileNames$sim[f])
# }
#
# # ===================================
# # Extraire les variables uniques
# # ===================================
# param1vars <- sort(unique(sapply(results, function(x) x$sim_effect$param_values[[param1Name]])))
# param2vars <- sort(unique(sapply(results, function(x) x$sim_effect$param_values[[param2Name]])))
# models <- unique(sapply(results, function(x) x$sim_effect$model))
# simus <- unique(sapply(results, function(x) x$replicat))
#
# # ===================================
# # Fonction d'extraction d'une simulation
# # ===================================
# select_results <- function(results, param1, param2, model, simu) {
#   for (i in seq_along(results)) {
#     r <- results[[i]]
#     if (r$sim_effect$param_values[[param1Name]] == param1 &&
#         r$sim_effect$param_values[[param2Name]] == param2 &&
#         r$sim_effect$model == model &&
#         r$replicat == simu) {
#
#       mediators <- r$sim_effect$mediators
#       res_feat <- r$res_effect$effects$univariate$ACME$feat
#       effet_A <- as.numeric(r$sim_effect$A_effect[mediators])
#       effet_B <- as.numeric(r$sim_effect$B_effect[mediators])
#       sign_simu <- sign(effet_A * effet_B)
#       sign_simu[sign_simu == -1] <- "-"
#       sign_simu[sign_simu == 1] <- "+"
#
#       res_med <- as.numeric(r$res_effect$effects$univariate$ACME$est)
#       sign_res <- sign(res_med)
#       sign_res[sign_res == -1] <- "-"
#       sign_res[sign_res == 1] <- "+"
#
#       names_all_probes <- union(names(mediators), res_feat)
#       data <- data.frame(
#         names_all_probes = names_all_probes,
#         simulation = rep("NS", length(names_all_probes)),
#         resultats = rep("NC", length(names_all_probes))
#       )
#
#       data$simulation[names_all_probes %in% names(mediators)] <- sign_simu
#       data$resultats[names_all_probes %in% res_feat] <- sign_res
#       return(data)
#     }
#   }
#   return(NULL)
# }
#
# # ===================================
# # Construction de toutes les combinaisons
# # ===================================
# all_results <- data.frame()
# for (param1var in param1vars) {
#   for (param2var in param2vars) {
#     for (model in models) {
#       for (s in simus) {
#         data <- select_results(results, param1var, param2var, model, s)
#         if (!is.null(data)) {
#           data$param1var <- param1var
#           data$param2var <- param2var
#           data$model <- model
#           all_results <- rbind(all_results, data)
#         }
#       }
#     }
#   }
# }
#
# # ===================================
# # Calcul des proportions dans la matrice de confusion
# # ===================================
# plot_data <- all_results %>%
#   group_by(simulation, resultats, param1var, param2var, model) %>%
#   summarise(Freq = n() / length(simus), .groups = "drop") %>%
#   mutate(Freq = ifelse(simulation == "NS" & resultats == "NC", NA, Freq))
#
# # ===================================
# # Affichage final
# # ===================================
# p <- ggplot(plot_data, aes(x = simulation, y = resultats, fill = Freq)) +
#   geom_tile(color = "white", na.rm = TRUE) +
#   geom_text(aes(label = sprintf("%.1f%%", 100 * Freq)), color = "black", size = 3, na.rm = TRUE) +
#   scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
#   facet_grid(param2var ~ param1var + model, labeller = label_both) +
#   theme_bw() +
#   coord_equal() +
#   theme(
#     strip.text = element_text(size = 8),
#     axis.text = element_text(size = 6),
#     axis.title = element_text(size = 6),
#     plot.title = element_text(size = 10)
#   ) +
#   labs(title = "Matrice de confusion normalisée", fill = "Proportion")
#
# # ===================================
# # Sauvegarde du graphique
# # ===================================
# ggsave(output_filename1, p, width = 14, height = 10)
#
# library(ggplot2)
# library(dplyr)
# library(stringr)
# library(tidyr)
#
# # ========================
# # PARAMÈTRES D'ANALYSE
# # ========================
#
# param1Name <- param1
# param2Name <- param2
#
# # ========================
# # LECTURE DU CSV D'ENTRÉE
# # ========================
# rdsFileNames <- read.csv(input_rds_files, header = FALSE, col.names = c("step2", "mediators"))
# rdsFileNames$sim <- as.numeric(str_match(rdsFileNames$step2, ".*_sim_([0-9]+).*")[, 2])
# rdsFileNames <- rdsFileNames[order(rdsFileNames$step2),]
#
# # ========================
# # CONSTRUCTION DU TABLEAU GLOBAL
# # ========================
# all_results <- data.frame()
#
# for (i in seq_len(nrow(rdsFileNames))) {
#   step2 <- readRDS(rdsFileNames$step2[i])
#   sim <- readRDS(rdsFileNames$mediators[i])
#   replicat <- rdsFileNames$sim[i]
#
#   p1val <- sim$param_values[[param1Name]]
#   p2val <- sim$param_values[[param2Name]]
#   model <- sim$model
#   method <- sim$method  # facultatif
#
#   # Effets simulés
#   causal_probes <- sim$mediators
#   effet_A <- as.numeric(sim$A_effect[causal_probes])
#   effet_B <- as.numeric(sim$B_effect[causal_probes])
#   sign_sim <- sign(effet_A * effet_B)
#   sign_sim[sign_sim == -1] <- "-"
#   sign_sim[sign_sim == 1] <- "+"
#
#   # Effets estimés
#   selected_probes <- step2$effects$univariate$ACME$feat
#   acme <- step2$effects$univariate$ACME$est
#   sign_est <- sign(as.numeric(acme))
#   sign_est[sign_est == -1] <- "-"
#   sign_est[sign_est == 1] <- "+"
#
#   # Union des sondes
#   all_probes <- union(names(sign_sim), selected_probes)
#
#   # ✅ Correction ici : utilisation de rep() pour toutes les colonnes
#   df <- data.frame(
#     probe_id = all_probes,
#     simulation = rep("NC", length(all_probes)),
#     resultats = rep("NS", length(all_probes)),
#     param1var = rep(p1val, length(all_probes)),
#     param2var = rep(p2val, length(all_probes)),
#     model = rep(model, length(all_probes)),
#     replicat = rep(replicat, length(all_probes))
#   )
#
#   # Ajout des vrais signes simulés et estimés
#   probes_in_sim <- df$probe_id %in% names(sign_sim)
#   df$simulation[probes_in_sim] <- sign_sim[df$probe_id[probes_in_sim]]
#
#   probes_in_est <- df$probe_id %in% selected_probes
#   df$resultats[probes_in_est] <- sign_est[df$probe_id[probes_in_est]]
#
#   # Ajout au tableau global
#   all_results <- bind_rows(all_results, df)
# }
#
# # ========================
# # AGGRÉGATION PAR COMBINAISON
# # ========================
# agg_results <- all_results %>%
#   group_by(simulation, resultats, param1var, param2var, model) %>%
#   summarise(Freq = n() / n_distinct(replicat), .groups = "drop") %>%
#   mutate(Freq = ifelse(simulation == "NC" & resultats == "NS", NA, Freq))
#
# # ========================
# # VISUALISATION
# # ========================
# p <- ggplot(agg_results, aes(x = simulation, y = resultats, fill = Freq)) +
#   geom_tile(color = "white", na.rm = TRUE) +
#   geom_text(aes(label = sprintf("%.1f%%", 100 * Freq)), color = "black", size = 3, na.rm = TRUE) +
#   scale_fill_gradient(low = "white", high = "steelblue", na.value = "white") +
#   facet_grid(model ~ param2var + param1var, labeller = label_both) +
#   coord_equal() +
#   theme_bw() +
#   labs(
#     title = "Matrice de confusion moyenne par combinaison de paramètres",
#     x = "Effets simulés",
#     y = "Effets estimés",
#     fill = "Proportion"
#   ) +
#   theme(
#     axis.text = element_text(size = 7),
#     strip.text = element_text(size = 8),
#     plot.title = element_text(size = 10)
#   )
#
# # ========================
# # EXPORT
# # ========================
# ggsave(output_filename1, p, width = 16, height = 10)
# p


# Juste pour un fichier au hasard

#
# i <- 1
# step2 <- readRDS(rdsFileNames$step2[i])
# sim <- readRDS(rdsFileNames$mediators[i])
#
# str(step2)  # important : contient-il $effects$univariate$ACME$feat et $est ?
# str(sim)    # important : contient-il $mediators, $A_effect, $B_effect, $param_values ?
#
# print(sim$model)
# print(sim$param_values)
# print(names(sim$mediators))  # sondes causales
# print(step2$effects$univariate$ACME$feat)  # sondes sélectionnées
#
#
#
# library(ggplot2)
#
# # ======= CONFIGURATION =======
# input_rds_files = "~/Datas/projects/thema_surv/simulations/v20/csv/n_nb_causal_probes_p_1000_n_n_pcp_nb_causal_probes_mA_2_mB_2_sA_0.1_sB_0.1_ro_0.5_aggregated_step2_results.csv"
# i <- 1  # index dans le CSV
# rdsFileNames <- read.csv(input_rds_files, header = FALSE, col.names = c("step2", "mediators"))
# step2 <- readRDS(rdsFileNames$step2[i])
# sim <- readRDS(rdsFileNames$mediators[i])
# replicat <- as.numeric(gsub(".*_sim_([0-9]+).*", "\\1", rdsFileNames$step2[i]))
#
# # ======= RÉCUPÉRATION DES DONNÉES =======
# causal_probes <- sim$mediators
# effet_A <- sim$A_effect[causal_probes]
# effet_B <- sim$B_effect[causal_probes]
# sign_sim <- sign(effet_A * effet_B)
# sign_sim[sign_sim == -1] <- "-"
# sign_sim[sign_sim == 1] <- "+"
#
# selected_probes <- step2$effects$univariate$ACME$feat
# acme <- step2$effects$univariate$ACME$est
# sign_est <- sign(acme)
# sign_est[sign_est == -1] <- "-"
# sign_est[sign_est == 1] <- "+"
#
# # ======= UNION DES SONS =======
# all_probes <- union(names(causal_probes), selected_probes)
# inter = intersect(names(causal_probes), selected_probes)
# df <- data.frame(
#   probe_id = all_probes,
#   simulation = rep("NC", length(all_probes)),
#   resultats = rep("NS", length(all_probes))
# )
#
# df$simulation[df$probe_id%in%names(causal_probes)] <- sign_sim
#
# df$resultats[df$probe_id %in% selected_probes] <- sign_est
#
# # ======= TABLE DE CONFUSION SIMPLE =======
# print(table(df$simulation, df$resultats))
#
# # ======= VISUALISATION SIMPLE =======
# confmat <- as.data.frame(table(df$simulation, df$resultats))
# colnames(confmat) <- c("simulation", "resultats", "Freq")
#
# ggplot(confmat, aes(x = simulation, y = resultats, fill = Freq)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = Freq), color = "black") +
#   theme_minimal() +
#   labs(title = paste("Matrice de confusion - replicat", replicat))
#
#

