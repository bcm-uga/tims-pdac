#################
### tims-pdac ###
#################

# Master script to run simulation analysis


##############
## Simulated data generation
##############
source("simulate_survival_mediation_data.R")
set.seed(s+42)

res = survival_data_simulation(n = n,
                               p = p,
                               nb.causal.probes = nb_causal_probes,
                               rho = rho,
                               sd.A = sd_A,
                               mean.A = mean_A,
                               sd.B = sd_B,
                               mean.B = mean_B,
                               causal_probes_overlap=causal_probes_overlap,
                               lambda_time = lambda_time)
param_values = list(n = n,
                    p = p,
                    nb_causal_probes = nb_causal_probes,
                    rho = rho,
                    sd_A = sd_A,
                    mean_A = mean_A,
                    sd_B = sd_B,
                    mean_B = mean_B,
                    causal_probes_overlap=causal_probes_overlap,
                    lambda_time = lambda_time)
loc_chr13 = readRDS("loc_chr13.rds")

loc_chr13_sort = sort(loc_chr13)[1:p]
loc_CpG = loc_chr13_sort[1:p]
meta_data = list(loc_CpG = loc_CpG)

generated_result = list(M=res$M, time=res$time, status=res$status, X=res$X, X_bin=res$X_bin, A_effect=res$A, B_effect=res$B, mediators=res$mediators, file=res$file, param_values=param_values, meta_data=meta_data)

saveRDS(result, filename)


##############
## HDMAX 2 Step 1 Parametric
##############

source("surv_parametric.R")

df = generated_result

param_values = df$param_values
param_values$model = "param"

#exposure = as.vector(df$X)
exposure = as.vector(df$X_bin)
survival_time = as.vector(df$time)
censoring_status = as.vector(as.numeric(df$status))
M = as.matrix(df$M)
K = 5

res = run_AS_surv_param(exposure,
                        survival_time,
                        censoring_status,
                        M, 
                        K,
                        covar = NULL,
                        suppl_covar = NULL,
                        each_var_pval = FALSE)


step1_result = list(res = res, param_values = param_values)





##############
## HDMAX 2 Step 1 Aalen
##############



source("surv_aalen.R")

df = generated_result

param_values = df$param_values
param_values$model = "aalen"

#exposure = as.vector(df$X)
exposure = as.vector(df$X_bin)
survival_time = as.vector(df$time)
censoring_status = as.vector(as.numeric(df$status))
M = as.matrix(df$M)
K = 5


res = run_AS_surv_aalen(exposure,
                        survival_time,
                        censoring_status,
                        M, 
                        K,
                        covar = NULL,
                        suppl_covar = NULL,
                        each_var_pval = FALSE)

step1_result = list(res = res, param_values = param_values)




##############
## HDMAX 2 Step 1 HIMA
##############

source("HIMA_last_version.R")

df = generated_result
param_values = df$param_values
param_values$model = "hima"

#exposure = as.vector(df$X)
exposure = as.vector(df$X_bin)
survival_time = as.vector(df$time)
censoring_status = as.vector(as.numeric(df$status))
M = as.matrix(df$M)
K = 5
#nb_select = unique(df$nb_causal_probes)

source("surv_parametric.R")

res_lfmm = lfmm2_med(input = M,
                     env = exposure, 
                     K= K, 
                     lambda = 1e-5,
                     effect.sizes = FALSE)


library(glmnet)
library(HDMT)
res_hima = hima_survival(exposure, M, survival_time, censoring_status, FDRcut = 1, topN = 100)

#transform res hima in hdmax2 output
res = list()
res$med = res_hima
res$max2_pvalues = res_hima$pmax
names(res$max2_pvalues) = res_hima$Index
res$AS_1$U = res_lfmm$U
res$input = list(
  exposure_input = exposure,
  expo_var_types = "double",
  expo_var_ids = "univariate",
  covar = NULL, 
  suppl_covar= NULL,
  survival_time_input = survival_time, 
  censoring_status_input = censoring_status
)


step1_result = list(res = res, param_values = param_values)



#############
## F1 score
#############
res_simu = generated_result
res_step1_content = step1_result
res_step1 = res_step1_content$res
param_values = res_step1_content$param_values


mediators = names(res_simu$mediators)

selected_med =  names(sort(res_step1$max2_pvalue)[1:(length(mediators))])


TP = length(intersect(mediators, selected_med)) 
# False Positive 
FP = length(selected_med) - TP

FN = length(mediators) - TP
precision = TP / (TP+FP)
recall = TP/(TP+FN)
F1_score_value = 2*(precision*recall)/(precision+recall)



res = data.frame(F1_score_value, TP, FN, FP, precision, recall, param_values)
saveRDS(res, output_filename)



##############
## HDMAX 2 Step 2 parametric
##############  
source("surv_parametric.R")

hdmax2_step1 = step1_result
hdmax2_step1 = hdmax2_step1$res
simu_var = readRDS(input_filename2)
M = as.matrix(simu_var$M)

m = M[, names(sort(hdmax2_step1$max2_pvalue))[1:simu_var$param_values$nb_causal_probes]]
step2_result = estimate_effect_surv_param(hdmax2_step1, m)



##############
## HDMAX 2 Step 2 Aalen
##############  

source("surv_aalen.R")

hdmax2_step1 = step1_result
hdmax2_step1 = hdmax2_step1$res
simu_var = readRDS(input_filename2)

M = as.matrix(simu_var$M)

## top 10 selection approach
m = M[, names(sort(hdmax2_step1$max2_pvalues))[1:simu_var$param_values$nb_causal_probes]]
step2_result = estimate_effect_surv_aalen(hdmax2_step1, m = m)

saveRDS(res, output_filename)

##############
## Step 2 HIMA
##############  

if (model == "hima") {
  hdmax2_step1 = step1_result
  hdmax2_step1 = hdmax2_step1$res
  step2_result = hdmax2_step1
  
}else {
  source("HIMA_last_version.R")
  library(glmnet)
  library(HDMT)
  hdmax2_step1 = step1_result
  hdmax2_step1 = hdmax2_step1$res
  simu_var = readRDS(input_filename2)
  M = as.matrix(simu_var$M)
  
  exposure = as.vector(simu_var$X_bin)
  survival_time = as.vector(simu_var$time)
  censoring_status = as.vector(as.numeric(simu_var$status))
  
  
  ##top 10 selection approach
  m = M[, names(sort(hdmax2_step1$max2_pvalues))[1:simu_var$param_values$nb_causal_probes]]
  m = as.matrix(m)
  step2_result = hima_survival(exposure, M, survival_time, censoring_status, FDRcut = 1, topN = length(simu_var$mediators))
  
}

###################################



