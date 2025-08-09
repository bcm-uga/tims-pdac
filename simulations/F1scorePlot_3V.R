
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

rdsFileNames = read.csv(input_rds_files, header = FALSE)


# Chargement et préparation des données
df = data.frame()
for (f in seq_along(rdsFileNames[, 1])) {
  f1score = readRDS(rdsFileNames[f,1])
  df = rbind(df, f1score)
}


# Garder les colonnes utiles
df = subset(df, select = c("F1_score_value", "TP", param1, param2, param3, "model"))
colnames(df) = c("F1Score", "TP", "param1", "param2", "param3", "model")

# Transformation en facteurs pour l'affichage
df$param1 <- factor(df$param1)
df$param2 <- factor(df$param2)
df$param3 <- factor(df$param3)
df$model  <- factor(df$model)




p <- ggplot(df, aes(x = model, y = F1Score, fill = factor(param3))) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.7)) +
  facet_grid(param2 ~ param1) +
  theme_bw() +
  labs(
    title = "Comparaison des F1-scores selon les modèles de survie",
    x = "Modèle de survie",
    y = "F1 score",
    fill = "Effectif"
  ) +
  ylim(0, 1) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 6),
    legend.position = "right"
  )

ggsave(output_filename1, p)
