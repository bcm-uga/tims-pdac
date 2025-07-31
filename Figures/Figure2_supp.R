################
### Figure 2 ###
################

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)
library(stringr)

param1 = "mean_A"
param2 = "mean_B"
param3 = "n"

# first panel
### Load and process data

rdsFileNames = read.csv("~pittionf/Datas/projects/thema_surv/simulations/v2-03/csv/mean_A_mean_B_n_p_10000_n_n_pcp_100_mA_mean_A_mB_mean_B_sA_0.1_sB_0.1_ro_0.5_overlap_0.8_lambda_0.1_aggregated_results.csv", header = FALSE)
df = data.frame()
for (f in seq_along(rdsFileNames[, 1])) {
  f1score = readRDS(rdsFileNames[f,1])
  df = rbind(df, f1score)
}


# keep usefull columns
df = subset(df, select = c("F1_score_value", "TP", param1, param2, param3, "model"))
colnames(df) = c("F1Score", "TP", "param1", "param2", "param3", "model")


df$param1 <- factor(df$param1)
df$param2 <- factor(df$param2)
df$param3 <- factor(df$param3)
df$model  <- factor(df$model)

# effects selection
keep1 <- df$param1 %in% c(0.1, 2)
keep2 <- df$param2 %in% c(0.01, 2)

df_sub <- df[keep1 & keep2, ]
df_sub <- droplevels(df_sub)
df_sub$model <- factor(
  df_sub$model,
  levels = c(
    "param","aalen","hima"
  )
)


# Création des fonctions de labelling
lab_param1 <- function(x) paste0("mean alpha = ", x)  
lab_param2 <- function(x) paste0("mean beta = ", x)  



big_theme <- theme(
  plot.title   = element_text(size = 16, face = "bold"),   
  axis.title   = element_text(size = 14),                  
  axis.text    = element_text(size = 12, angle = 45, hjust = 1),                  
  strip.text   = element_text(size = 13, face = "bold"),  
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 12),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(),
  panel.grid.minor   = element_blank()
)

p1 <- ggplot(df_sub, aes(
  x    = model,
  y    = TP,
  fill = factor(param3)
)) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.6,
    position      = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(color = factor(param3)),
    position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2),
    size       = 1,
    alpha      = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8, name = "Samples") +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  facet_grid(
    param2 ~ param1,
    labeller = labeller(
      param1 = lab_param1,
      param2 = lab_param2
    )
  ) +
  theme_bw() +
  labs(
    title = "True Positive Selection Across Mediation Approaches",
    x     = "Survival mediation methods",
    y     = "True Positive"
  ) +
  ylim(0, 100) +
  big_theme

#print(p1)


ggsave("figures/step1_results_supp.pdf", p1, width = 10.5, height = 7)



# Second panel
### Load and process data

rdsFileNames = read.csv("~pittionf/Datas/projects/thema_surv/simulations/v2-03/csv/mean_A_mean_B_n_p_10000_n_n_pcp_100_mA_mean_A_mB_mean_B_sA_0.1_sB_0.1_ro_0.5_overlap_0.8_lambda_0.1_aggregated_step2_results.csv", header = FALSE, col.names = c("step2", "mediators"))

sim = as.numeric(str_match(rdsFileNames$step2, ".*_sim_([012457]+).*")[,2])
rdsFileNames = cbind(rdsFileNames, sim)
rdsFileNames <- rdsFileNames[order(rdsFileNames$step2),]

param1 = "mean_A"
param2 = "mean_B"
param3 = "n"

results = vector(mode = "list")
for (f in 1:length(rdsFileNames$step2)) {
  res_effect = readRDS(rdsFileNames$step2[f])
  #res_effects[[f]] = res_effect
  sim_effect = readRDS(rdsFileNames$mediators[f])
  #sim_effects[[f]] = sim_effect
  results[[f]] = list(res_effect = res_effect ,sim_effect=sim_effect, replicat = sim[f])
  #results$sim_effects[[f]] = sim_effect
}

param1Name = param1
param2Name = param2
param3Name = param3


param1vars = param2vars = param3vars = models = methods = model_methods = simus = c()


for (r in 1:length(results)){
  param1var = results[[r]]$sim_effect$param_values[[param1Name]]
  param2var = results[[r]]$sim_effect$param_values[[param2Name]]
  param3var = results[[r]]$sim_effect$param_values[[param3Name]]
  param1vars = c(param1vars, param1var)
  param2vars = c(param2vars, param2var)
  param3vars = c(param3vars, param3var)
  model = results[[r]]$sim_effect$model
  method = results[[r]]$sim_effect$method
  model_method = paste0(model,"/", method)
  simu = results[[r]]$replicat
  models = c(models, model)
  methods = c(methods, method)
  simus = c(simus, simu)
  model_methods = c(model_methods, model_method)
}

param1vars = unique(param1vars)
param1vars = sort(param1vars)
param2vars = unique(param2vars)
param2vars = sort(param2vars)
param3vars = unique(param3vars)
param3vars = sort(param3vars)
models = unique(models)
simus = unique(simus)
methods = unique(methods)
model_methods = unique(model_methods)


model_methods_df = data.frame(model =c("aalen","aalen", "aalen","hima", "hima","hima", "param", "param", "param"), method = c("hdmax2a","hdmax2p", "hima","hdmax2a","hdmax2p", "hima","hdmax2a","hdmax2p", "hima"), model_method = model_methods)

## function for plot df compute good precision score

select_results <- function(results, param1Value, param2Value, param3Value, model, method, simu, alpha = 0.05) {
  for (i in seq_along(results)) {
    r <- results[[i]]
    
    p1_match <- r$sim_effect$param_values[[param1Name]] == param1Value
    p2_match <- r$sim_effect$param_values[[param2Name]] == param2Value
    p3_match <- r$sim_effect$param_values[[param3Name]] == param3Value
    model_match <- r$sim_effect$model == model
    simu_match <- r$replicat == simu
    method_match =r$sim_effect$method == method
    
    if (p1_match && p2_match&& p3_match && model_match && simu_match && method_match) {
      
      mediators_name <- names(r$sim_effect$mediators)
      
      effet_A <- as.numeric(r$sim_effect$A_effect[r$sim_effect$mediators])
      effet_B <- as.numeric(r$sim_effect$B_effect[r$sim_effect$mediators])
      sign_simu <- sign(effet_A * effet_B)
      sign_simu[sign_simu == -1] <- "-"
      sign_simu[sign_simu == 1] <- "+"
      
      if(method == "hima" && model == "hima"){
        res_med <- as.numeric(r$res_effect$med$IDE)
        res_feat <- r$res_effect$med$Index
        sign_res <- sign(res_med)
        sign_res[sign_res == -1] <- "-"
        sign_res[sign_res == 1] <- "+"
        significante = r$res_effect$med$pmax  <= alpha/r$sim_effect$param_values$nb_causal_probes
      } else if(method == "hima" && model != "hima"){
        res_med <- as.numeric(r$res_effect$IDE)
        res_feat <- r$res_effect$Index
        sign_res <- sign(res_med)
        sign_res[sign_res == -1] <- "-"
        sign_res[sign_res == 1] <- "+"
        significante = r$res_effect$pmax  <= alpha/r$sim_effect$param_values$nb_causal_probes
      } else{
        res_med <- as.numeric(r$res_effect$effects$univariate$ACME$est)
        res_feat <- r$res_effect$effects$univariate$ACME$feat
        sign_res <- sign(res_med)
        sign_res[sign_res == -1] <- "-"
        sign_res[sign_res == 1] <- "+"
        significante = r$res_effect$effects$univariate$ACME$pval <= alpha/r$sim_effect$param_values$nb_causal_probes
      }
      
      all_probes <- union(mediators_name, res_feat)
      
      data <- data.frame(
        names_all_probes = all_probes,
        simulation = rep("NC", length(all_probes)),
        resultats = rep("NS", length(all_probes)),
        significante = rep(FALSE, length(all_probes))
      )
      
      data$simulation[all_probes %in% mediators_name] <- sign_simu
      data$resultats[all_probes %in% res_feat] <- sign_res
      data$significante[all_probes %in% res_feat] = significante
      return(data)
    }
  }
  
}

## plot df

scores            <- c()
model_methods     <- c()
param1_cumul      <- c()
param2_cumul      <- c()
param3_cumul      <- c()


param1vars = c(0.1, 2)
param2vars = c(0.01, 2)

# param1vars = c(0, 0.1, 2)
# param2vars = c(0, 0.01, 2)


for (i in seq_len(nrow(model_methods_df))) {
  model        <- model_methods_df$model[i]
  method       <- model_methods_df$method[i]
  model_method <- model_methods_df$model_method[i]
  
  for (param1 in param1vars) {
    for (param2 in param2vars) {
      for (param3 in param3vars) {
        for (s in seq_along(simus)) {
          
          # extraction des résultats filtrés
          data_plot <- select_results(
            results     = results,
            param1Value = param1,
            param2Value = param2,
            param3Value = param3,
            model       = model,
            method      = method,
            simu        = simus[s],
            alpha       = 0.05
          )
          
          # calcul du score
          score <- 0
          for (j in seq_len(nrow(data_plot))) {
            sim  <- data_plot$simulation[j]
            res  <- data_plot$resultats[j]
            sig  <- data_plot$significante[j]
            
            if (method %in% c("hdmax2a", "hima")) {
              if (sim == res && sig) score <- score + 1
            } else {
              if (sim != res && sig) score <- score + 1
            }
          }
          
          scores        <- c(scores,        score / nrow(data_plot))
          model_methods <- c(model_methods, model_method)
          param1_cumul  <- c(param1_cumul,  param1)
          param2_cumul  <- c(param2_cumul,  param2)
          param3_cumul  <- c(param3_cumul,  param3)
        }
      }
    }
  }
}




df <- data.frame(
  model_method = model_methods,
  param1       = param1_cumul,
  param2       = param2_cumul,
  param3       = param3_cumul,
  score        = scores
)
df$model_method <- factor(
  df$model_method,
  levels = c(
    "param/hdmax2p","aalen/hdmax2p","hima/hdmax2p","param/hdmax2a","aalen/hdmax2a","hima/hdmax2a","param/hima",
    "aalen/hima","hima/hima"
  )
)



## plot

theme_set(theme_bw(base_size = 12))   # <- police
big_theme <- theme(
  plot.title   = element_text(size = 16, face = "bold"),   
  axis.title   = element_text(size = 14),                  
  axis.text    = element_text(size = 12, angle = 45, hjust = 1),                 
  strip.text   = element_text(size = 13, face = "bold"),  
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 12),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(),
  panel.grid.minor   = element_blank()
)


col_labeller  <- function(x) paste0("Mean alpha = ", x)  # colonnes (param1)
row_labeller  <- function(x) paste0("Mean beta = ", x)  # lignes   (param2)


dodge_w <- 0.8

p2 <- ggplot(df,
             aes(x = model_method,
                 y = score,
                 fill = factor(param3))) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.6,
    position      = position_dodge2(width = dodge_w, padding = 0.1),
    alpha         = 0.75
  ) +
  geom_point(
    aes(color = factor(param3)),
    position = position_jitterdodge(
      dodge.width  = dodge_w,
      jitter.width = 0.15,
      jitter.height = 0
    ),
    size  = 0.5,
    alpha = 0.55
  ) +
  facet_grid(
    rows     = vars(param2),
    cols     = vars(param1),
    labeller = labeller(.rows = row_labeller,
                        .cols = col_labeller)   # ← les bons labellers
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8, name = "Samples") +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  
  scale_x_discrete(
    breaks = c("param/hdmax2p","aalen/hdmax2p","hima/hdmax2p",
               "param/hdmax2a","aalen/hdmax2a","hima/hdmax2a",
               "param/hima", "aalen/hima","hima/hima"),
    labels = c("param_s1/param_s2","aalen_s1/param_s2","hima/param_s2",
               "param_s1/aalen_s2","aalen_s1/aalen_s2","hima/aalen_s2",
               "param_s1/hima","aalen_s1/hima","hima/hima")
  ) +
  theme_bw() +
  labs(
    title = "Good precision score distribution by A and B mean effect",
    x     = "Method",
    y     = "Good precision score",
    fill  = "Samples",
    color = "Samples"
  ) +
  ylim(0, 1) +
  big_theme 

#print(p2)

ggsave("figures/step2_res_plot_supp.pdf",p2 , width = 10.5, height = 7)


