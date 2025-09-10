library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)
param1 = "mean_A"
param2 = "mean_B"
param3 = "n"




###########################################
# TP among second reg pval2 with 0 B effect
###########################################


# load data
simus = read.csv("simulation/v2-10/csv/all_simu.csv", header = FALSE)
step1s = read.csv("simulation/v2-10/csv/all_step1s.csv", header = FALSE)


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



lab_param1 <- function(x) paste0("mean alpha = ", x)
lab_param2 <- function(x) paste0("mean beta = ", x)

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
      meanA = lab_param1,
      meanB = lab_param2
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

ggsave("figures/step1_results_top_pval2_0_Beffect.pdf", p_pval2_0B, width = 10.5, height = 5)