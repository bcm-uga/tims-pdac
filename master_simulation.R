############################################################
### tims-pdac — DEMO (non-executable)
### Goal: illustrate the workflow for a set of parameters
############################################################

## ===============================
## 0) Preamble (sources & parameters)
## ===============================
# We assume that the functions are defined in the sourced scripts.
# This is only to illustrate the logic, not to run the analysis.

# -- Data generation
source("simulations/simulate_survival_mediation_data.R")

# -- Survival models (parametric & Aalen)
source("simulations/surv_parametric.R")
source("simulations/surv_aalen.R")

# -- HIMA method (selection + harmonized output)
source("simulations/HIMA_last_version.R")

# Parameters here are placeholders (not defined),
# just to keep the workflow readable:
# n, p, nb_causal_probes, rho, sd_A, mean_A, sd_B, mean_B,
# causal_probes_overlap, lambda_time, and seed "s".
set.seed(s + 42)


## =================================================
## 1) Generate simulated survival mediation data
## =================================================
# Inputs -> simulator ; Outputs -> list() "generated_result"
res <- survival_data_simulation(
  n = n,
  p = p,
  nb.causal.probes = nb_causal_probes,
  rho = rho,
  sd.A = sd_A,
  mean.A = mean_A,
  sd.B = sd_B,
  mean.B = mean_B,
  causal_probes_overlap = causal_probes_overlap,
  lambda_time = lambda_time
)

# Example metadata (CpG coordinates, truncated to p)
loc_chr13      <- readRDS("loc_chr13.rds")
loc_chr13_sort <- sort(loc_chr13)[1:p]
loc_CpG        <- loc_chr13_sort[1:p]

param_values <- list(
  n = n,
  p = p,
  nb_causal_probes = nb_causal_probes,
  rho = rho,
  sd_A = sd_A,
  mean_A = mean_A,
  sd_B = sd_B,
  mean_B = mean_B,
  causal_probes_overlap = causal_probes_overlap,
  lambda_time = lambda_time
)

meta_data <- list(loc_CpG = loc_CpG)

generated_result <- list(
  M          = res$M,                  # mediator matrix
  time       = res$time,               # survival time
  status     = res$status,             # censoring indicator
  X          = res$X,                  # continuous exposure (if available)
  X_bin      = res$X_bin,              # binary exposure (used below)
  A_effect   = res$A,                  # simulated A effects
  B_effect   = res$B,                  # simulated B effects
  mediators  = res$mediators,          # ground truth causal mediators
  file       = res$file,               # file name (if saving)
  param_values = param_values,
  meta_data    = meta_data
)


## ==========================================
## 2) HDMAX2 — Step 1 (mediator screening)
##    -> Three variants: Parametric / Aalen / HIMA
## ==========================================

# ------- Common inputs for Step 1 -------
df                <- generated_result
param_values_step <- df$param_values
exposure          <- as.vector(df$X_bin)      # binary exposure (choice here)
survival_time     <- as.vector(df$time)
censoring_status  <- as.vector(as.numeric(df$status))
M                 <- as.matrix(df$M)
K                 <- 5                         # number of latent factors (example)
covar <- NULL; suppl_covar <- NULL             # no covariates in this demo

# 2A) Step 1 — Parametric
param_values_param       <- param_values_step; param_values_param$model <- "param"
res_param <- run_AS_surv_param(
  exposure, survival_time, censoring_status, M, K,
  covar = covar, suppl_covar = suppl_covar, each_var_pval = FALSE
)
step1_param <- list(res = res_param, param_values = param_values_param)

# 2B) Step 1 — Aalen 
param_values_aalen       <- param_values_step; param_values_aalen$model <- "aalen"
res_aalen <- run_AS_surv_aalen(
  exposure, survival_time, censoring_status, M, K,
  covar = covar, suppl_covar = suppl_covar, each_var_pval = FALSE
)
step1_aalen <- list(res = res_aalen, param_values = param_values_aalen)

# 2C) Step 1 — HIMA (harmonized output to HDMAX2 format)
param_values_hima        <- param_values_step; param_values_hima$model <- "hima"

# (Optional) latent factors via lfmm2 to align outputs
res_lfmm <- lfmm2_med(input = M, env = exposure, K = K, lambda = 1e-5, effect.sizes = FALSE)

# HIMA selection (survival)
res_hima_raw <- hima_survival(
  exposure, M, survival_time, censoring_status,
  FDRcut = 1, topN = 100
)

# Wrap into HDMAX2-like format
res_hima <- list()
res_hima$med            <- res_hima_raw
res_hima$max2_pvalues   <- res_hima_raw$pmax
names(res_hima$max2_pvalues) <- res_hima_raw$Index
res_hima$AS_1$U         <- res_lfmm$U
res_hima$input <- list(
  exposure_input           = exposure,
  expo_var_types           = "double",
  expo_var_ids             = "univariate",
  covar                    = covar,
  suppl_covar              = suppl_covar,
  survival_time_input      = survival_time,
  censoring_status_input   = censoring_status
)
step1_hima <- list(res = res_hima, param_values = param_values_hima)


## ==========================================
## 3) Intermediate evaluation — F1 score
##    (Compare selection vs ground truth)
## ==========================================
# Example: compute F1 for one Step 1 result (parametric).
# Same logic applies for aalen/hima.

res_simu            <- generated_result
res_step1_content   <- step1_param                 # <<< choose param / aalen / hima
res_step1           <- res_step1_content$res
param_values_f1     <- res_step1_content$param_values

mediators_true <- names(res_simu$mediators)        # ground truth causal mediators
# Selection: sort by $max2_pvalues, keep as many as true mediators
selected_med <- names(sort(res_step1$max2_pvalues))[1:length(mediators_true)]

TP <- length(intersect(mediators_true, selected_med))
FP <- length(selected_med) - TP
FN <- length(mediators_true) - TP

precision <- TP / (TP + FP)
recall    <- TP / (TP + FN)
F1_score  <- 2 * (precision * recall) / (precision + recall)

f1_table <- data.frame(
  F1_score = F1_score, TP = TP, FN = FN, FP = FP,
  precision = precision, recall = recall,
  model = param_values_f1$model
)


## ==========================================
## 4) HDMAX2 — Step 2 (effect estimation)
##    -> Based on mediators selected in Step 1
## ==========================================

# Take one Step 1 result:
hdmax2_step1_out <- step1_param$res   # <<< choose param / aalen / hima

# Mediator matrix from simulated variables
simu_var <- generated_result
M2       <- as.matrix(simu_var$M)

# Selection strategy for Step 2:
nb_keep <- simu_var$param_values$nb_causal_probes
m_sel   <- M2[, names(sort(hdmax2_step1_out$max2_pvalues))[1:nb_keep]]

# ===== 4A) Step 2 — Parametric =====
step2_param <- estimate_effect_surv_param(hdmax2_step1_out, m = m_sel)

# ===== 4B) Step 2 — Aalen =====
step2_aalen <- estimate_effect_surv_aalen(hdmax2_step1_out, m = m_sel)

# ===== 4C) Step 2 — HIMA =====
if (param_values_f1$model == "hima") {
  step2_hima <- hdmax2_step1_out
} else {
  step2_hima <- hima_survival(exposure, M2, survival_time, censoring_status, FDRcut = 1, topN   = length(simu_var$mediators))
}


