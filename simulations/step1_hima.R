#load last hima version last release not effective
source("HIMA_last_version.R")

df = readRDS(input_filename)
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


result = list(res = res, param_values = param_values)
saveRDS(result, output_filename)


