library(cmdstanr)

## Load data
y = as.matrix(read.csv("./Data/Ukraine_y.csv", header = F))
x_cov = as.matrix(read.csv("./Data/Ukraine_x_cov.csv", header = F))
x_cov = x_cov[,-c(1, 9)] ## Remove intercept and last column
age = x_cov[,2]
age = scale(age)
age_sq = age^2
age_sq = age_sq - mean(age_sq)
x_cov = x_cov[,-2]
x_cov = scale(x_cov, center = TRUE, scale = FALSE)

z_cov = as.matrix(read.csv("./Data/resplevel_Ukraine_z_cov.csv", header = F))
z_cov = scale(z_cov, center = TRUE, scale = FALSE)



N = 46510000
known <- c(4088438, 935153, 1328606, 3966956, 889443, 2993846, 1869098, 258619, 150989, 144130,104186, 47587, 278195, 222884, 762877,
           323863, 108511, 178364, 273200, 472657, 69471, 487148)
names(known)<- c("M20-30", "M15-17", "M70+", "F20-30", "F15-17", "F70+", "Kids", "Moldovans", "Romanians", "Poles", "Jews","Roma",
                 "Invalids", "MedicalDoctors", "Died2007", "Pavlo", "Prisoners2007", "DivorcedMen2007","Policemen", "Birth2007", "PhD", "Nurses")
n.known = length(known)

#########
## For removing 9 bad subpopulations
keep.known = c(1, 3, 4, 6, 7, 13, 15, 16, 17, 18, 20)
keep.vec = c(keep.known, c(23, 24, 25, 26))
y = y[,keep.vec]
z_cov = z_cov[,keep.vec]
known = known[keep.known]
n.known = length(known)
##
##########
N.i = nrow(y)
N.k = ncol(y)

stan.data = list(
  n_i = nrow(y),
  n_k = ncol(y),
  x_size = ncol(x_cov),
  x_cov = x_cov,
  z_cov = z_cov,
  age = c(age),
  age_sq = c(age_sq),
  y = y
)

stan.model = cmdstan_model("./Stan_Files/Correlated_Stan.stan")

model.fit = stan.model$sample(data = stan.data,
                              chains = 3,
                              parallel_chains = 3,
                              iter_warmup = 1500,
                              iter_sampling = 500,
                              refresh = 100,
                              init = 0)


logdi.est = model.fit$draws(variables = "logdi",
                           format = "draws_matrix")
sigma.di.est = model.fit$draws(variables = "sigma_di",
                            format = "draws_matrix")
prev.est = model.fit$draws(variables = "prevalence",
                                    format = "draws_matrix")
beta.est = model.fit$draws(variables = "beta_par",
                           format = "draws_matrix")
alpha.est = model.fit$draws(variables = "alpha_par",
                           format = "draws_matrix")
age.est = model.fit$draws(variables = "age_par",
                            format = "draws_matrix")
age.sq.est = model.fit$draws(variables = "age_sq_par",
                          format = "draws_matrix")
tauN = model.fit$draws(variables = "tau_N",
                             format = "draws_matrix")
bias = model.fit$draws(variables = "bias",
                       format = "draws_matrix")
Corr = model.fit$draws(variables = "Corr",
                           format = "draws_matrix")
prev.mean = model.fit$draws(variables = "prev_mean",
                       format = "draws_matrix")


saveRDS(logdi.est, file = "./Results/Ukraine_correlated_logdi.rds")
saveRDS(sigma.di.est, file = "./Results/Ukraine_correlated_sigma_di.rds")
saveRDS(prev.est, file = "./Results/Ukraine_correlated_prevalence.rds")
saveRDS(beta.est, file = "./Results/Ukraine_correlated_beta_par.rds")
saveRDS(alpha.est, file = "./Results/Ukraine_correlated_alpha_par.rds")
saveRDS(age.est, file = "./Results/Ukraine_correlated_age_par.rds")
saveRDS(age.sq.est, file = "./Results/Ukraine_correlated_age_sq_par.rds")
saveRDS(tauN, file = "./Results/Ukraine_correlated_tau_N.rds")
saveRDS(bias, file = "./Results/Ukraine_correlated_bias.rds")
saveRDS(Corr, file = "./Results/Ukraine_correlated_Corr.rds")
saveRDS(prev.mean, file = "./Results/Ukraine_correlated_prev_mean.rds")

saveRDS(model.fit, file = "./Results/Ukraine_correlated_fit.rds")
## Not all of these are used, but may be useful for further analysis