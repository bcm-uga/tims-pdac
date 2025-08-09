source("surv_aalen.R")

df = readRDS(input_filename)

# ####################"
# df = result
# ################
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

result = list(res = res, param_values = param_values)

saveRDS(result, output_filename)


