args = commandArgs(trailingOnly=TRUE)
model = args[1]
method = args[2]
input = args[3]
output = args[4]

simu = readRDS(input)
data = list(
  mediators = simu$mediators,
  A_effect = simu$A_effect,
  B_effect = simu$B_effect,
  param_values = simu$param_values,
  model = model,
  method = method
)
saveRDS(data, output)
