
res_simu = readRDS(input_filename)
res_step1_content = readRDS(input_filename2)
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


