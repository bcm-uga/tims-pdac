med_tobacco_dnam_CpG_param = readRDS("results/02_tcga_med_tobacco_dnam_V2_K8_corrected.rds")
med_tobacco_dnam_CpG_hima = readRDS("results/02_tcga_med_tobacco_dnam_hima.rds")

vec1 = names(sort(med_tobacco_dnam_CpG_param$hdmax2_step1_param$max2_pvalues)[1:100])

pval = as.matrix(med_tobacco_dnam_CpG_hima$hdmax2_step1_hima$pvalue_SIS)
pval = t(pval)
cpg = med_tobacco_dnam_CpG_hima$hdmax2_step1_hima$out_result$Index
rownames(pval) = cpg
vec2 = names(sort(med_tobacco_dnam_CpG_hima$hdmax2_step1_hima$pvalue_SIS)[1:100])
length(intersect(vec1, cpg))
