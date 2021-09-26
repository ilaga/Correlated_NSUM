
library(cmdstanr)

# Fit it on Ukraine -----------------------------------------------

y = as.matrix(read.csv("./Data/Ukraine_y.csv", header = F))
N = 46510000
known <- c(4088438, 935153, 1328606, 3966956, 889443, 2993846, 1869098, 258619, 150989, 144130,104186, 47587, 278195, 222884, 762877,
           323863, 108511, 178364, 273200, 472657, 69471, 487148)
names(known)<- c("M20-30", "M15-17", "M70+", "F20-30", "F15-17", "F70+", "Kids",
                 "Moldovans", "Romanians", "Poles", "Jews","Roma",
                 "Invalids", "MedicalDoctors", "Died2007", "Pavlo",
                 "Prisoners2007", "DivorcedMen2007","Policemen", "Birth2007",
                 "PhD", "Nurses")
n.known = length(known)

## Removing 9 bad subpopulations
keep.known = c(1, 3, 4, 6, 7, 13, 15, 16, 17, 18, 20)
keep.vec = c(keep.known, c(23, 24, 25, 26))
## Order of unknown is FSW, MSW, MSM, and IDU
y = y[,keep.vec]
known = known[keep.known]
n.known = length(known)

N.i = nrow(y)
N.k = ncol(y)

stan.model = cmdstan_model("./Stan_Files/Overdispersed_Stan.stan")

stan.data = list(
  n_i = N.i,
  n_k = N.k,
  y = y
)

overdispersed.Ukraine = stan.model$sample(data = stan.data,
                              chains = 3,
                              parallel_chains = 3,
                              iter_warmup = 1500,
                              iter_sampling = 500,
                              refresh = 100,
                              init = 0)

betas.est = overdispersed.Ukraine$draws(variables = "betas",
                            format = "draws_matrix")
alphas.est = overdispersed.Ukraine$draws(variables = "alphas",
                                format = "draws_matrix")
omegas.est = 1/ overdispersed.Ukraine$draws(variables = "inv_omegas",
                                 format = "draws_matrix") ## Want omegas, so take inverse

saveRDS(betas.est, file = "./Results/Ukraine_overdispersed_betas.rds")
saveRDS(alphas.est, file = "./Results/Ukraine_overdispersed_alphas.rds")
saveRDS(omegas.est, file = "./Results/Ukraine_overdispersed_omegas.rds")
saveRDS(overdispersed.Ukraine, file = "./Results/Ukraine_overdispersed_fit.rds")
## Not used elsewhere, but may be useful for further analysis
