
library(ggplot2)
library(stringr)
library(patchwork)
library(grid)
library(gridExtra)


# recup des res step2 pour unr methode faire les compte dans les ^
#input_rds_files = "/home/pittionf/Datas/projects/thema_surv/simulations/v23/csv/n_nb_causal_probes_causal_probes_overlap_p_1000_n_n_pcp_nb_causal_probes_mA_2_mB_2_sA_0.1_sB_0.1_ro_0.5_overlap_causal_probes_overlap_aggregated_step2_results.csv"
#input_rds_files = "/home/pittionf/Datas/projects/thema_surv_archive/v4/csv/nb_causal_probes_mean_A_p_1000_n_500_pcp_nb_causal_probes_mA_mean_A_mB_2_sA_0.5_sB_0.5_ro_0.5_aggregated_step2_results.csv"
#input_rds_files = "/home/pittionf/Datas/projects/thema_surv/simulations/v13/csv/mean_A_mean_B_n_p_500_n_n_pcp_40_mA_mean_A_mB_mean_B_sA_0.1_sB_0.1_ro_0.5_aggregated_step2_results.csv"

rdsFileNames = read.csv(input_rds_files, header = FALSE, col.names = c("step2", "mediators"))
sim = as.numeric(str_match(rdsFileNames$step2, ".*_sim_([0-9]+).*")[,2])
rdsFileNames = cbind(rdsFileNames, sim)
rdsFileNames <- rdsFileNames[order(rdsFileNames$step2),]


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
param3Name = param3

param1vars = param2vars = param3vars = models = simus = c()

 for (r in seq_along(results)){
   simu = results[[r]]$replicat
   param1var = results[[r]]$sim_effect$param_values[[param1Name]]
   param2var = results[[r]]$sim_effect$param_values[[param2Name]]
   param3var = results[[r]]$sim_effect$param_values[[param3Name]]
   param1vars = c(param1vars, param1var)
   param2vars = c(param2vars, param2var)
   param3vars = c(param3vars, param3var)
   model = results[[r]]$sim_effect$model
   models = c(models, model)
   simus = c(simus, simu)
 }
   
param1vars = unique(param1vars)
param1vars = sort(param1vars)
param2vars = unique(param2vars)
param2vars = sort(param2vars)
param3vars = unique(param3vars)
param3vars = sort(param3vars)
models = unique(models)
simus = unique(simus)

# select_results = function(results, param1 , param2, param3, model, simu){
#   for(i in 1:length(results)) {
#     if(results[[i]]$sim_effect$param_values[[param1Name]]==param1 && results[[i]]$sim_effect$param_values[[param2Name]]==param2 && results[[i]]$sim_effect$param_values[[param3Name]]==param3 && results[[i]]$sim_effect$model==model && results[[i]]$replicat==simu) {
#       mediators = results[[i]]$sim_effect$mediators
#       res_feat = results[[i]]$res_effect$effects$univariate$ACME$feat
#       inter = intersect(names(mediators), res_feat)
#       sign_simu = sign(results[[i]]$sim_effect$A_effect[results[[i]]$sim_effect$mediators]*results[[i]]$sim_effect$B_effect[results[[i]]$sim_effect$mediators])
#       sign_simu[sign_simu == -1] = "-"
#       sign_simu[sign_simu == 1] = "+"
#       sign_res = sign(results[[i]]$res_effect$effects$univariate$ACME$est)
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
#     }
#   }
# }
# 

select_results <- function(results, param1Value, param2Value, param3Value, model, simu) {
  for (i in seq_along(results)) {
    r <- results[[i]]

    p1_match <- r$sim_effect$param_values[[param1Name]] == param1Value
    p2_match <- r$sim_effect$param_values[[param2Name]] == param2Value
    p3_match <- r$sim_effect$param_values[[param3Name]] == param3Value
    model_match <- r$sim_effect$model == model
    simu_match <- r$replicat == simu

    if (p1_match && p2_match&& p3_match && model_match && simu_match) {
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
  stop(sprintf("Aucun match trouvé pour param1 = %s, param2 = %s,param3 = %s, model = %s, simu = %s", param1Value, param2Value, param3Value, model, simu))

}

# Initialisation de la liste des graphiques
plots <- list()
param3var = param3vars[1]
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
          param3Value = param3var,
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


# Initialisation de la liste des graphiques
plots <- list()
param3var = param3vars[2]
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
          param3Value = param3var,
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
ggsave(output_filename2, final_plot, width = 10, height = 8)


#
# # Initialisation de la liste des graphiques
# plots <- list()
# param3var = param3vars[3]
# model_colors <- c("param" = "lightblue", "aalen" = "pink")
#
#
#
# # Stockage dans une liste nommée : [[param2]][[param1]][[model]]
# for (param2var in param2vars) {
#   for (param1var in param1vars) {
#     for (model in models) {
#
#       data_sign <- data.frame()
#
#       for (s in seq_along(simus)) {
#         data <- select_results(
#           results = results,
#           param1Value = param1var,
#           param2Value = param2var,
#           param3Value = param3var,
#           model = model,
#           simu = simus[s]
#         )
#         data_sign <- rbind(data_sign, data)
#       }
#
#       confusion_matrix <- table(data_sign$simulation, data_sign$resultats)
#       plot_data <- as.data.frame(confusion_matrix)
#       colnames(plot_data) <- c("simulation", "resultats", "Freq")
#       plot_data$Freq <- plot_data$Freq / length(simus)
#       plot_data$Freq[plot_data$simulation == "NS" & plot_data$resultats == "NC"] <- NA
#
#       p <- ggplot(plot_data, aes(x = simulation, y = resultats, fill = Freq)) +
#         geom_tile(color = "white", na.rm = TRUE) +
#         scale_fill_gradient(
#           high = scales::muted(model_colors[model]),
#           low = model_colors[model],
#           na.value = "white"
#         ) +
#         coord_equal() +
#         geom_text(aes(label = round(Freq, 0)), size = 4, na.rm = TRUE) +
#         theme_minimal(base_size = 8) +
#         labs(
#           x = "Simulation",
#           y = "Résultats",
#           title = paste("Model:", model),
#           subtitle = paste(param1, "=", param1var, "|", param2, "=", param2var)
#         ) +
#         theme(
#           plot.title = element_text(face = "bold", size = 9),
#           plot.subtitle = element_text(size = 8),
#           axis.text = element_text(size = 10, face = "bold")
#         )
#
#       plots[[paste(param2var, param1var, model, sep = "_")]] <- p
#     }
#   }
# }
#
# # Réorganiser les plots pour les empiler par lignes/colonnes
# # => chaque ligne = un param2, chaque colonne = un param1
# grid_plots <- list()
# for (param2val in param2vars) {
#   row_plots <- list()
#   for (param1val in param1vars) {
#     row_plots <- c(row_plots, list(
#       plots[[paste(param2val, param1val, "param", sep = "_")]],
#       plots[[paste(param2val, param1val, "aalen", sep = "_")]]
#     ))
#   }
#   grid_plots <- c(grid_plots, row_plots)
# }
#
# # Créer la grille 2x2 ou plus selon le nombre de paramètres
# final_plot <- arrangeGrob(
#   grobs = grid_plots,
#   nrow = length(param2vars),
#   ncol = 2 * length(param1vars)  # 2 colonnes par valeur de param1 (param + aalen)
# )
#
# # Affichage et export
# grid.newpage()
# grid.draw(final_plot)
# ggsave(output_filename3, final_plot, width = 10, height = 8)
#
#