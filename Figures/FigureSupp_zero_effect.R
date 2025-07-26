library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)
param1 = "mean_A"
param2 = "mean_B"
param3 = "n"

##########
# 0 B effect plot
#########
rdsFileNames = read.csv("~pittionf/Datas/projects/thema_surv/simulations/v2-10/csv/mean_A_mean_B_n_p_10000_n_n_pcp_100_mA_mean_A_mB_mean_B_sA_0.1_sB_0_ro_0_overlap_0.8_lambda_0.1_aggregated_results.csv", header = FALSE)

df = data.frame()
for (f in seq_along(rdsFileNames[, 1])) {
  f1score = readRDS(rdsFileNames[f,1])
  df = rbind(df, f1score)
}

df = subset(df, select = c("F1_score_value", "TP", param1, param2, param3, "model"))
colnames(df) = c("F1Score", "TP", "param1", "param2", "param3", "model")

df$param1 <- factor(df$param1)
df$param2 <- factor(df$param2)
df$param3 <- factor(df$param3)
df$model  <- factor(df$model)

df_sub =df

df_sub$model <- factor(
  df_sub$model,
  levels = c(
    "param","aalen","hima"
  )
)

lab_param1 <- function(x) paste0("mean A = ", x)
lab_param2 <- function(x) paste0("mean B = ", x)

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

p_zeroB <- ggplot(df_sub, aes(
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

#print(p_zeroB)

ggsave("figures/step1_results_0_Beffect.pdf", p_zeroB, width = 10.5, height = 7)

##########""
# 0 A effect plot
#########
rdsFileNames = read.csv("~pittionf/Datas/projects/thema_surv/simulations/v2-11/csv/mean_A_mean_B_n_p_10000_n_n_pcp_100_mA_mean_A_mB_mean_B_sA_0_sB_0.1_ro_0.5_overlap_0.8_lambda_0.1_aggregated_results.csv", header = FALSE)
df = data.frame()
for (f in seq_along(rdsFileNames[, 1])) {
  f1score = readRDS(rdsFileNames[f,1])
  df = rbind(df, f1score)
}

df = subset(df, select = c("F1_score_value", "TP", param1, param2, param3, "model"))
colnames(df) = c("F1Score", "TP", "param1", "param2", "param3", "model")

df$param1 <- factor(df$param1)
df$param2 <- factor(df$param2)
df$param3 <- factor(df$param3)
df$model  <- factor(df$model)

df_sub =df

df_sub$model <- factor(
  df_sub$model,
  levels = c(
    "param","aalen","hima"
  )
)
library(ggplot2)
library(viridis)

lab_param1 <- function(x) paste0("mean A = ", x)
lab_param2 <- function(x) paste0("mean B = ", x)

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

p_zeroA <- ggplot(df_sub, aes(
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
    param1 ~ param2,
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

print(p_zeroA)

ggsave("figures/step1_results_0_Aeffect.pdf", p_zeroA, width = 10.5, height = 7)


###########################################
# TP among second reg pval2 with 0 B effect
###########################################


# load data
simus = read.csv("~/Datas/projects/thema_surv/simulations/v2-10/csv/all_simu.csv", header = FALSE)
step1s = read.csv("~/Datas/projects/thema_surv/simulations/v2-10/csv/all_step1s.csv", header = FALSE)


# function
parseSimuVars <- function(file) {
  n=str_match(file, "_n_([0-9]+)_.*")[,2]
  mA=str_match(file, "_mA_([0-9.]+)_.*")[,2]
  mB=str_match(file, "_mB_([0-9.]+)_.*")[,2]
  sim=str_match(file, "_sim_([0-9]+)*")[,2]
  return (list(n=n,mA=mA, mB=mB,sim=sim))
}

computeRes <- function(simufile, step1file) {
  simu = readRDS(simufile)
  step1 = readRDS(step1file)
  pval2 = step1$res$AS_2$pval
  idxmed = names(simu$mediators)
  idxb = names(sort(pval2))[1:100]
  return( length(intersect(idxmed, idxb)))
}

# df for plot
c = 1
df=data.frame(model = character(), n = character(), meanA = character() , meanB = character(), sim = character(), res = numeric() )
for (s in 1:30) {
  simufile  = simus[[1]][s]
  aalenfile = step1s[[1]][c]
  paramfile = step1s[[1]][c+2]
  simuvars = parseSimuVars(simufile)
  resaalen = computeRes(simufile , aalenfile)
  resparam = computeRes(simufile , paramfile)
  new_line = data.frame(model="aalen", n=simuvars$n, meanA = simuvars$mA , meanB = simuvars$mB, sim = simuvars$sim , res = resaalen)
  df = rbind(df, new_line)
  new_line = data.frame(model="param", n=simuvars$n, meanA = simuvars$mA , meanB = simuvars$mB, sim = simuvars$sim , res = resparam)
  df = rbind(df, new_line)
  c = c+3
  
}



lab_param1 <- function(x) paste0("mean A = ", x)
lab_param2 <- function(x) paste0("mean B = ", x)

big_theme <- theme(
  plot.title   = element_text(size = 16, face = "bold"),   
  axis.title   = element_text(size = 14),                  
  axis.text    = element_text(size = 14),                  
  strip.text   = element_text(size = 13, face = "bold"),  
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 12),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(),
  panel.grid.minor   = element_blank()
)
p_pval2_0B <- ggplot(df, aes(
  x    = model,
  y    = res,
  fill = factor(n)
)) +
  geom_boxplot(
    outlier.shape = NA,
    width         = 0.6,
    position      = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(color = factor(n)),
    position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2),
    size       = 1,
    alpha      = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8, name = "Samples") +
  scale_color_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  facet_grid(
    meanB ~ meanA,
    labeller = labeller(
      param1 = lab_param1,
      param2 = lab_param2
    )
  ) +
  theme_bw() +
  labs(
    title = "True Positive among second mediation regression",
    x     = "Survival mediation methods",
    y     = "True Positive"
  ) +
  ylim(0, 100) +
  big_theme

print(p_pval2_0B)


ggsave("figures/step1_results_top_pval2_0_Beffect.pdf", p_pval2_0B, width = 10.5, height = 5)
