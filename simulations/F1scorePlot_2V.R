library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

rdsFileNames = read.csv(input_rds_files, header = FALSE)

df = data.frame()
 for (f in 1:length(rdsFileNames[,1])) {
    f1score = readRDS(rdsFileNames[f,1])
    df=rbind(df, f1score)
 }


#corriger method vs model apres avoir retourné step1
#df = subset(df, select = c("F1_score_value", param1, param2,"method"))# "model")
df = subset(df, select = c("F1_score_value", param1, param2, "model"))
colnames(df) = c("F1Score", "param1", "param2","model")


# Convertir les colonnes en facteurs pour ggplot
df$model <- as.factor(df$model)
df$param2 <- as.factor(df$param2)

# Création du scatter plot

data <- df %>%
  group_by(param2, model) %>%
  mutate(id = row_number()) %>%
  ungroup()

# Appliquer pivot_wider pour avoir une colonne pour chaque méthode
data_wide <- data %>%
  pivot_wider(names_from = model, values_from = F1Score)

p = ggplot(data_wide, aes(x = param1, y = param2)) +
  geom_point(aes(x = aalen, y = param)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_grid(param2 ~ param1) +
  xlim(0, 1)+
  ylim(0, 1)+
  theme_bw() +
  labs(x = "F1Score, survival model = aalen", y = "F1Score, survival model =  param", color = "nb of samples", shape = "nb of samples") +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 6))

ggsave(output_filename1, p)



