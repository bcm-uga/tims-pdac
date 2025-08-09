source("new_test_simu.R")
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

result = list(M=res$M, time=res$time, status=res$status, X=res$X, X_bin=res$X_bin, A_effect=res$A, B_effect=res$B, mediators=res$mediators, file=res$file, param_values=param_values, meta_data=meta_data)

saveRDS(result, filename)