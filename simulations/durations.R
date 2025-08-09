library(stringr)
library(ggplot2)

#inputfilename = "v2-02/summary_execution_times.csv"

durations = read.csv(inputfilename)
durations = cbind(durations, p=str_match(durations$step, "_p_([0-9]+)_.*")[,2])
durations = cbind(durations, n=str_match(durations$step, "_n_([0-9]+)_.*")[,2])
durations = cbind(durations, pcp=str_match(durations$step, "_pcp_([0-9]+)_.*")[,2])
durations = cbind(durations, mA=str_match(durations$step, "_mA_([0-9.]+)_.*")[,2])
durations = cbind(durations, mB=str_match(durations$step, "_mB_([0-9.]+)_.*")[,2])
durations = cbind(durations, sA=str_match(durations$step, "_sA_([0-9.]+)_.*")[,2])
durations = cbind(durations, sB=str_match(durations$step, "_sB_([0-9.]+)_.*")[,2])
durations = cbind(durations, ro=str_match(durations$step, "_ro_([0-9.]+)_.*")[,2])
durations = cbind(durations, overlap=str_match(durations$step, "_overlap_([0-9.]+)_.*")[,2])
durations = cbind(durations, lambda=str_match(durations$step, "_lambda_([0-9.]+)_.*")[,2])
durations = cbind(durations, model=str_match(durations$step, "_model_([a-z]+)_.*")[,2])

if (hasName(durations, "total")) {
  durations = cbind(durations, time=durations$total)
} else {
  durations = cbind(durations, time=durations$step1)
}

durations$log_time = log10(durations$time)

model_colors <- c("param" = "blue", "aalen" = "red", "hima" = "green")

plot = ggplot(durations, aes(x =  n , y = log_time, color = model)) +
  geom_point(cex = 0.3) +
  facet_grid(mA~mB)+
  
  scale_color_manual(               
    values = model_colors,              
    name   = "Modèle",             
    breaks = names(model_colors)        
  )+
  theme_bw()
ggsave(outputfileGrid, plot)

library(ggplot2)



plot2 <- ggplot(durations,
                aes(x = factor(n), y = log_time,
                    fill  = model,          # couleur des violons
                    color = model)) +       # même couleur pour les points
  ## 1. Violons
  geom_violin(
    position = position_dodge(width = 0.8),
    trim     = TRUE,# modif
    width    = 0.9,
    alpha    = 0.7
  ) +
  ## 2. Points individuels (jitter)
  geom_point(
    position = position_jitterdodge(
      dodge.width  = 0.8,   # même décalage que les violons
      jitter.width = 0.15
    ),
    size  = 1.3,
    alpha = 0.7
  ) +
  labs(
    title = "Execution Time Distributions (log-scale)",
    x     = "Sample Size (n)",
    y     = "log10(Execution Time)",
    fill  = "Method",
    color = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(size = 14, face = "bold")
  )+ expand_limits(y = 0) # modif

#print(plot2)

ggsave(outputfileGlobal, plot2, width = 7, height = 5)


#ggsave("figures/durations_global.pdf", plot2, width = 7, height = 5)
