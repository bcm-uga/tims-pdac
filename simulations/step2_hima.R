if (model == "hima") {
  hdmax2_step1 = readRDS(input_filename)
  hdmax2_step1 = hdmax2_step1$res
  res = hdmax2_step1
  saveRDS(res, output_filename)
}else {
  source("HIMA_last_version.R")
  library(glmnet)
  library(HDMT)
  hdmax2_step1 = readRDS(input_filename)
  hdmax2_step1 = hdmax2_step1$res
  simu_var = readRDS(input_filename2)
  M = as.matrix(simu_var$M)

  exposure = as.vector(simu_var$X_bin)
  survival_time = as.vector(simu_var$time)
  censoring_status = as.vector(as.numeric(simu_var$status))


  ##top 10 selection approach
  m = M[, names(sort(hdmax2_step1$max2_pvalues))[1:simu_var$param_values$nb_causal_probes]]
  m = as.matrix(m)
  res = hima_survival(exposure, M, survival_time, censoring_status, FDRcut = 1, topN = length(simu_var$mediators))
  saveRDS(res, output_filename)
}

