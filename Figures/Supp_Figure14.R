######################
### Supp Figure 14 ###
######################
library(ggplot2)

### Load and process data

# Load the TCGA-PAAD data
tcga_data = readRDS("tcga_pdac_mediation/results/01_tcga_data_expo_deconv.rds") # Load the TCGA-PAAD data

# Load the immune cell-type estimation
datIMM = read.csv2("tcga_pdac_mediation/results/03_tcga_consensus_deconv_immune_cells.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datIMM[] <- lapply(datIMM, function(x) as.numeric(as.character(x)))
datIMM$all = rowSums(datIMM)
colnames(datIMM) = c("Macrophages", "B cells", "T cells", "NK", "DCs", "Tot. Imm." )

# Load the significant AMR results
tobacco_AMR <- readRDS("tcga_pdac_mediation/results/03_tcga_significative_from_top50_tobacco_AMR_fdr0_05_V2_K8_corrected.rds")
AMR_names = tobacco_AMR$AMR_info$AMR

# Load the AMR results (methylation value per AMR)
datAMR = read.csv2("tcga_pdac_mediation/results/03_tcga_AMR_mean_meth_top50_fdr0_05_V2_K8_corrected.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",") 
datAMR[] <- lapply(datAMR, function(x) as.numeric(as.character(x)))

# Process the smoking variable from TCGA-PAAD metadata
smoking = ifelse(tcga_data$tobacco==0, 'Non-smoker', 'Smoker')
names(smoking) = rownames(tcga_data$M)
labels = smoking[rownames(datAMR)]
smoking_status = as.numeric(as.factor(labels))

#Latent factor composition and correlation
res_med = readRDS("tcga_pdac_mediation/results/02_tcga_med_tobacco_dnam_V2_K8_corrected.rds")
LFs = res_med$hdmax2_step1_param$AS_1$U
colnames(LFs) = paste0("LF_", c("A", "B", "C", "D", "E", "F", "G", "H"))
pairs_list = readRDS("tcga_pdac_mediation/results/05_signLFs_by_pairs-A-I.rds")

### Run analysis

##### For Immune infiltration

# Compute Pearson correlation matrix between LFs and datIMM
cor_matrix_LF <- cor(LFs, datIMM, use = "pairwise.complete.obs", method = "pearson")

# Function to get p-value from Pearson correlation test between two vectors
get_pval <- function(x, y) cor.test(x, y, use = "pairwise.complete.obs", method = "pearson")$p.value

# Vectorized computation of p-values for all combinations of columns in LFs and datIMM
pval_matrix_LF <- outer(
  1:ncol(LFs), 1:ncol(datIMM),
  Vectorize(function(i, j) get_pval(LFs[, i], datIMM[, j]))
)

# Set row and column names of the p-value matrix to match LFs and datIMM column names
dimnames(pval_matrix_LF) <- list(colnames(LFs), colnames(datIMM))

# Bonferroni adjusted significance threshold
alpha_adj <- 0.01 / (50+ 6) # 50 LFs + 6 immune variables

# Convert correlation matrix to a data frame for ggplot
df_plot <- as.data.frame(as.table(cor_matrix_LF))
colnames(df_plot) <- c("LF", "Immune", "Correlation")

# Add corresponding p-values to the data frame
df_plot$p.value <- mapply(function(lf, im) pval_matrix_LF[lf, im], df_plot$LF, df_plot$Immune)

# Mark correlations as strong if absolute correlation is greater than 0.4
df_plot$Strong <- abs(df_plot$Correlation) > 0.4

df_plot$significance <- ifelse(df_plot$p.value <= alpha_adj, "Significant", "Not Significant")

# Create dot plot with correlation magnitude as size, correlation value as color,

p_imm = ggplot(df_plot, aes(x = Immune, y = LF)) +
  geom_point(aes(size = abs(Correlation),
                 color = Correlation,
                 shape = significance),
             stroke = 1.2, fill = "white") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  scale_size(range = c(0, 5), name = "|Correlation|") +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 17)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlations between LFs and immune variables",
       shape = "pval < 0.01",
       x = "Immune",
       y = "LF")

##### For AMRs

# Compute Pearson correlation matrix between LFs and datAMR
cor_matrix_AMR <- cor(LFs, datAMR, use = "pairwise.complete.obs", method = "pearson")

# Function to get p-value from Pearson correlation test between two vectors
get_pval <- function(x, y) cor.test(x, y, use = "pairwise.complete.obs", method = "pearson")$p.value

# Vectorized computation of p-values for all combinations of columns in LFs and datAMR
pval_matrix_AMR <- outer(
  1:ncol(LFs), 1:ncol(datAMR),
  Vectorize(function(i, j) get_pval(LFs[, i], datAMR[, j]))
)

# Set row and column names of the p-value matrix to match LFs and datAMR column names
dimnames(pval_matrix_AMR) <- list(colnames(LFs), colnames(datAMR))

# Bonferroni adjusted significance threshold
alpha_adj <- 0.01 / (50+ 6) # 50 LFs + 6 immune variables

# Convert correlation matrix to a data frame for ggplot
df_plot_AMR <- as.data.frame(as.table(cor_matrix_AMR))
colnames(df_plot_AMR) <- c("LF", "AMR", "Correlation")

# Add corresponding p-values to the data frame
df_plot_AMR$p.value <- mapply(function(lf, amr) pval_matrix_AMR[lf, amr], df_plot_AMR$LF, df_plot_AMR$AMR)

# Mark correlations as strong if absolute correlation is greater than 0.4
df_plot_AMR$Strong <- abs(df_plot_AMR$Correlation) > 0.4

df_plot_AMR$significance <- ifelse(df_plot_AMR$p.value <= alpha_adj, "Significant", "Not Significant")


p_amr = ggplot(df_plot_AMR, aes(x = AMR, y = LF)) +
  geom_point(aes(size = abs(Correlation),
                 color = Correlation,
                 shape = significance),
             stroke = 1.2, fill = "white") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  scale_size(range = c(0, 5), name = "|Correlation|") +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 17)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Correlations between LFs and immune variables",
       shape = "pval < 0.01",
       x = "AMR",
       y = "LF")

### Combine and save the plot

combined_plot <- p_imm + p_amr + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("figures/supp_fig14.pdf", combined_plot,
       width = 14, height = 6, units = "in")
