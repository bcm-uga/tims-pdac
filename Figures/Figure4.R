################
### Figure 4 ###
################

library(survival)
library(survminer)
library(gridExtra)
library(ggplot2)  
library(ComplexHeatmap)

### Load and process data

# Load the TCGA-PAAD data
tcga_data = readRDS("real_data/results/01_tcga_data_expo_deconv.rds") # Load the TCGA-PAAD data

# Load the immune cell-type estimation
datIMM = read.csv2("real_data/results/03_tcga_consensus_deconv_immune_cells.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",")
datIMM[] <- lapply(datIMM, function(x) as.numeric(as.character(x)))
datIMM$all = rowSums(datIMM)
colnames(datIMM) = c("Macrophages", "B cells", "T cells", "NK", "DCs", "Tot. Imm." )

# Load the significant AMR results
tobacco_AMR <- readRDS("real_data/results/03_tcga_significative_from_top50_tobacco_AMR_fdr0_05_V2_K8_corrected.rds")
AMR_names = tobacco_AMR$AMR_info$AMR

# Load the AMR results (methylation value per AMR)
datAMR = read.csv2("real_data/results/03_tcga_AMR_mean_meth_top50_fdr0_05_V2_K8_corrected.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ",") 
datAMR[] <- lapply(datAMR, function(x) as.numeric(as.character(x)))

# Process the smoking variable from TCGA-PAAD metadata
smoking = ifelse(tcga_data$tobacco==0, 'Non-smoker', 'Smoker')
names(smoking) = rownames(tcga_data$M)
labels = smoking[rownames(datAMR)]
smoking_status = as.numeric(as.factor(labels))

#Latent factor composition and correlation
res_med = readRDS("real_data/results/02_tcga_med_tobacco_dnam_V2_K8_corrected.rds")
LFs = res_med$hdmax2_step1_param$AS_1$U
colnames(LFs) = paste0("LF_", c("A", "B", "C", "D", "E", "F", "G", "H"))
pairs_list = readRDS("real_data/results/05_signLFs_by_pairs-A-I.rds")


### PANEL A: TCGA-PAAD heatmap of deconvolution results

X_deconv = as.data.frame(t(tcga_data$X_deconv$Consensus))
X_deconv$TotalImmune = c(tcga_data$prop_immune)
range(X_deconv)
X_deconv = X_deconv[,c("TotalImmune","B cells","T cells","DCs","NK","Neutrophils","Macrophages",
                       "Endothelial","Fibroblasts",
                       "Cancer basal","Cancer classical")]
colnames(X_deconv) = c("Tot. Imm.","B cells","T cells","DCs","NK","Neutrophils","Macrophages",
                       "Endothelial","Fibroblasts",
                       "Cancer basal","Cancer classical")
					   
Z_score = apply(X_deconv, 2, function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T))
range(Z_score, na.rm=T)

# Annotation pre processing
gender = ifelse(tcga_data$gender == 1, "male", "female")
grade = tcga_data$grade
stopifnot(length(gender) == nrow(Z_score))
stopifnot(length(grade) == nrow(Z_score))

# Row annotation
row_anno = rowAnnotation(
  Gender = gender,
  Grade = grade,
  col = list(
    Gender = c("male" = "brown", "female" = "yellow"),
    Grade = structure(RColorBrewer::brewer.pal(length(unique(grade)), "Set2"),
                      names = unique(grade))
  ),
  show_annotation_name = TRUE
)


pdf("figures/fig4_panelA.pdf", width = 4, height =4)
ComplexHeatmap::Heatmap(Z_score[,colnames(Z_score) != 'Neutrophils'],
                        cluster_rows = T,
                        show_row_names = F,
                        cluster_columns = F,
						 left_annotation = row_anno)
dev.off()

### PANEL B: survival curves for immune cell types


# Run cox model and select significant immmune features

pval_thres = 0.15

res = apply(datIMM, 2, function(x) {
  model = survival::coxph(survival::Surv(tcga_data$time, tcga_data$status) ~ x + tcga_data$age + tcga_data$gender +  tcga_data$grade)
  summary(model)$coefficients[, "Pr(>|z|)"]
})

print(res)
immune_names = names(which(res[1,] <= pval_thres)) # select imm -> survival in the graph
length(immune_names)

# Plot Kaplan-Meier

selected_vars <- colnames(datIMM)[1:6]  
plots <- list()

time_in_months <- tcga_data$time / 30.44  # 1 month ≈ 30.44 jours

for (var in selected_vars) {
  
  # Dichotomize based on the median
  group <- ifelse(datIMM[[var]] >= median(datIMM[[var]], na.rm = TRUE), "High", "Low")
  group <- factor(group, levels = c("Low", "High"))
  plot_data <- data.frame(time = time_in_months, event = tcga_data$status, group = group)
  
  # Créer l'objet de survie
  surv_obj <- survival::Surv(time = time_in_months, event = tcga_data$status)
  
  # Kaplan-Meier survival estimation
  fit <- survival::survfit(surv_obj ~ group)
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = plot_data,
    risk.table = FALSE,
    pval = TRUE,
    title = var,
    legend.title = var,
    legend.labs = c("Low", "High"),
    palette = c("#E7B800", "#2E9FDF"),
    xlab = "Time (months)"
  )
  
  plots[[var]] <- p
}

#  Display the 6 survival curves in a 2x3 grid
combined_plot = arrange_ggsurvplots(plots, ncol = 3, nrow = 2)

ggsave("figures/fig4_panelB.pdf", plot = combined_plot, width = 8, height = 5, dpi = 300)

### PANEL C: Causal discovery

pval_thres = 0.1

# details of each model for publication:

prop100= datIMM[,"Tot. Imm."]*100 #to interpret the HR has increase of 1 unit (1% of immune infiltration) lead to HR likely to dire than persons with 1% less
mod = survival::coxph(survival::Surv(tcga_data$time, tcga_data$status) ~ prop100 + tcga_data$age + tcga_data$gender +  tcga_data$grade)
summary(mod)

prop100= datIMM[,"B cells"]*100 #to interpret the HR has increase of 1 unit (1% of immune infiltration) lead to HR likely to dire than persons with 1% less
mod = survival::coxph(survival::Surv(tcga_data$time, tcga_data$status) ~ prop100 + tcga_data$age + tcga_data$gender +  tcga_data$grade)
summary(mod)

prop100= datIMM[,"DCs"]*100 #to interpret the HR has increase of 1 unit (1% of immune infiltration) lead to HR likely to dire than persons with 1% less
mod = survival::coxph(survival::Surv(tcga_data$time, tcga_data$status) ~ prop100 + tcga_data$age + tcga_data$gender +  tcga_data$grade)
summary(mod)

# Latent factor assessment

surv_LF = apply(LFs, 2, function(x) {
  survival::coxph(survival::Surv(tcga_data$time, tcga_data$status) ~ x)
})
surv_LF

tob_LF = apply(LFs, 2, function(x) {
  cor.test(x, smoking_status)$p.value
})
tob_LF


# Step 1: Unconditional independence resting

## Tobacco-Survival association: Cox proportional hazards model

SURV_TOB = survival::coxph(survival::Surv(tcga_data$time, tcga_data$status) ~  smoking_status + tcga_data$age + tcga_data$gender +  tcga_data$grade )
summary(SURV_TOB)

## Tobacco-Immune associations: Linear regression models

IMMtot_TOB = lm(datIMM[ ,"Tot. Imm."]~smoking_status + tcga_data$age + tcga_data$gender +  tcga_data$grade )    
summary(IMMtot_TOB) 

IMMbcell_TOB = lm(datIMM[ ,"B cells"]~smoking_status + tcga_data$age + tcga_data$gender +  tcga_data$grade )    
summary(IMMbcell_TOB) 

IMMDCs_TOB = lm(datIMM[ ,"DCs"]~smoking_status + tcga_data$age + tcga_data$gender +  tcga_data$grade)
summary(IMMDCs_TOB) 


AMR_names = tobacco_AMR$AMR_info$AMR

num_res = list()

for (AMR in AMR_names){ 
  
  for (imm in immune_names) {
    
    
    pair = paste(AMR, "-", imm, sep = "")
    print(paste0("Testing AMR:", AMR, " and Immune:", imm))
    
    # Get significant LFs 
    LFs_names = pairs_list[[pair]]
    covar =   data.frame(age = tcga_data$age,
                         gender = tcga_data$gender,
                         grade =   tcga_data$grade
    )
    if (is.na(LFs_names[1])) {
      print(paste0("No significant LFs for pair: ", pair))
    } else {
      print(paste0("Significant LFs for pair: ", pair, " are: ", paste(LFs_names, collapse = ", ")))
      covar = cbind(covar, LFs[, LFs_names, drop = FALSE])
    }
    
    
    # Check consisty of T -> AMR -> S path
    
    #  Check the link between tobacco and the AMR 
    
    df_model <- data.frame(
      y = datAMR[, AMR],
      smoking_status = smoking_status,
      covar  
    )
    mod_tob <- lm(y ~ ., data = df_model)
    tob_sign =  summary(mod_tob)$coefficients[2,4]  < pval_thres 
    
    #  check the link between the AMR and survival
    
    df_model <- data.frame(
      y = datAMR[, AMR],
      time = tcga_data$time,
      status = tcga_data$status,
      covar)
    
    mod_AMR <- survival::coxph(
      survival::Surv(time, status) ~ y +  ., 
      data = df_model  
    )
    
    AMR_sign =  summary(mod_AMR)$coefficients[1, "Pr(>|z|)"] < pval_thres 
    
    # Step 2: Conditional independence testing
    
    if (tob_sign == TRUE & AMR_sign == TRUE) {
      
      print(paste0(" AMR:", AMR, " and Immune:", imm, " is kept for doanwstream analysis. "))
      # (Model 1) Which node coefficient will lose its significance between AMR and imm ?
      
      df_model <- data.frame(
        time = tcga_data$time,
        status = tcga_data$status,
        imm = datIMM[[imm]],
        amr = datAMR[[AMR]],
        covar
      )
      mod_part <- survival::coxph(
        survival::Surv(time, status) ~ imm + amr +  .,
        data = df_model
      )
      part_surv_sign = summary(mod_part)$coefficients[1:2, "Pr(>|z|)"] 
      
      
      # (Model 2) check the link between the immune and smoking status when controlling for AMR
      
      
      df_model <- data.frame(
        smoking_status = as.numeric(as.factor(labels)) -1,
        AMR = datAMR[,AMR],
        imm = datIMM[ ,imm],
        covar
      )
      mod = glm(smoking_status~ ., data = df_model, family = binomial(link = "logit"))
      summary(mod)
      part_imm_sign = summary(mod)$coefficients[2:3,4] # check if imm is significant
      
      num_res[[pair]] = c( part_surv_sign,
                           part_imm_sign)
      names(num_res[[pair]]) = c("part_IMM_sign_to_surv", 
                                 "part_AMR_sign_to_surv_when_IMM_controlled",
                                 "part_AMR_sign_to_smoking",
                                 "part_IMM_sign_to_smoking_when_AMR_controlled")
    } else {
      print(paste0(" AMR:", AMR, " and Immune:", imm, " is not kept for doanwstream analysis. "))
      #num_res[[pair]] = c(NA, NA, NA, NA)
      #names(num_res[[pair]]) = c("part_IMM_sign_to_surv", 
      #                            "part_AMR_sign_to_surv_when_IMM_controlled",
      #                           "part_AMR_sign_to_smoking",
      #                           "part_IMM_sign_to_smoking_when_AMR_controlled")
    }
    
  }
  
  
}
num_res = do.call(rbind, num_res)
num_res=data.frame(num_res)

# Directed acyclic graph 1
res_CAT1 = num_res[  num_res$part_IMM_sign_to_surv < pval_thres &
                       num_res$part_AMR_sign_to_surv_when_IMM_controlled < pval_thres, ]
res_CAT1
rownames(res_CAT1)

# Directed acyclic graph 2
res_CAT2 = num_res[  num_res$part_IMM_sign_to_surv < pval_thres &
                       num_res$part_AMR_sign_to_surv_when_IMM_controlled >= pval_thres, ]
rownames(res_CAT2)



plot_part_cor = function(mat){
  
  colnames(mat) = c("Model 1: S~I+A+C, (pval I)", 
                    "Model 1: S~I+A+C, (pval A)",
                    "Model 2: T~I+A+C, (pval A)",
                    "Model 2: T~I+A+C, (pval I)")
  
  df <- data.frame(
    cell_type = rep(rownames(mat), times = ncol(mat)),
    variable = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(as.matrix(mat))
  )
  
  df$abs_value <- abs(df$value)
  df$significance <- ifelse(df$value <= 0.1, "Significant", "Not Significant")
  
  
  p = ggplot(df, aes(x = variable, y = cell_type,
                     size = abs_value,
                     color = value,
                     shape = significance)) +
    geom_point() +
    scale_size(range = c(2, 8)) +
    scale_color_gradient(low = "darkblue", high = "lightblue") +
    scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 17)) +  # 16 = rond, 17 = triangle
    theme_minimal() +
    labs(
      title = "Partial correlation",
      x = "Linear models",
      y = "Pairs",
      color = "P-value",
      size = "Intensity",
      shape = "Significativity"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


pdf("figures/fig4_panelC1.pdf", width = 5, height =5)
plot_part_cor(mat = res_CAT1)
dev.off()

  pdf("figures/fig4_panelC2.pdf", width = 5, height =5)
plot_part_cor(mat = res_CAT2)
dev.off()
