source("surv_parametric.R")

hdmax2_step1 = readRDS(input_filename)
hdmax2_step1 = hdmax2_step1$res
simu_var = readRDS(input_filename2)
M = as.matrix(simu_var$M)

m = M[, names(sort(hdmax2_step1$max2_pvalue))[1:simu_var$param_values$nb_causal_probes]]
res = estimate_effect_surv_param(hdmax2_step1, m)
saveRDS(res, output_filename)

