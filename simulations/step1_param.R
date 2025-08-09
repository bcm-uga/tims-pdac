
source("surv_parametric.R")

df = readRDS(input_filename)

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

 
result = list(res = res, param_values = param_values)
saveRDS(result, output_filename)

