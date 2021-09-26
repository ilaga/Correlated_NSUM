
library(mvtnorm)
library(cmdstanr)

ind = ## Fill in simulation number
size.in = ## Fill in sample size
model = ## Fill in the model here from the following list
# c("Corr_x_z", "Corr_nox_z", "Corr_x_noz", "Uncorr_x_z", "Uncorr_nox_noz")


## Load data
y = as.matrix(read.csv("./Data/Ukraine_y.csv", header = F))
x_cov = as.matrix(read.csv("./Data/Ukraine_x_cov.csv", header = F))
x_cov = x_cov[,-c(1, 9)] ## Remove intercept and last column
age = x_cov[,2]
age = scale(age)
age_sq = age^2
age_sq = age_sq - mean(age_sq)
x_cov = x_cov[,-2]
z_cov = as.matrix(read.csv("./Data/resplevel_Ukraine_Z_cov.csv", header = F))
## Center z_cov by 3
z_cov = z_cov - 3

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


## Sample covariates
sample.ind = sample(1:N.i, size.in, replace = TRUE)
N.i = size.in
age = age[sample.ind]
age_sq = age_sq[sample.ind]
x_cov = x_cov[sample.ind,]
z_cov = z_cov[sample.ind,]

## Center
age = scale(age)
age_sq = age^2
age_sq = age_sq - mean(age_sq)
x_cov = scale(x_cov, center = TRUE, scale = FALSE)
z_cov = scale(z_cov, center = TRUE, scale = FALSE)



corr.prevalence = colMeans(readRDS("./Results/Ukraine_correlated_prevalence.rds"))
corr.corr.in = readRDS("./Results/Ukraine_correlated_Corr.rds")
corr.corr = array(NA, dim = c(nrow(corr.corr.in), N.k, N.k))
for(i in 1:nrow(corr.corr.in)){
  corr.corr[i,,] = matrix(corr.corr.in[i,], nrow = N.k, ncol = N.k)
}
corr.corr.ave = apply(corr.corr, c(2,3), mean)
corr.tau = colMeans(readRDS("./Results/Ukraine_correlated_tau_N.rds"))
corr.alpha = colMeans(readRDS("./Results/Ukraine_correlated_alpha_par.rds"))
corr.beta = colMeans(readRDS("./Results/Ukraine_correlated_beta_par.rds"))
corr.age = colMeans(readRDS("./Results/Ukraine_correlated_age_par.rds"))
corr.agesq = colMeans(readRDS("./Results/Ukraine_correlated_age_sq_par.rds"))
corr.sigma.di = mean(readRDS("./Results/Ukraine_correlated_sigma_di.rds"))


logdegree.sim = rnorm(N.i) * corr.sigma.di
mvt.mean = log(1 / sqrt(1 + corr.tau^2))
mvt.tau = sqrt(log(1 + corr.tau^2))
bias.sim = rmvnorm(N.i, mvt.mean, diag(mvt.tau) %*% corr.corr.ave %*% diag(mvt.tau))
y.sim = matrix(NA, nrow = N.i, ncol = N.k)
for(i in 1:N.i){
  y.sim[i,] = rpois(N.k, exp(logdegree.sim[i] + corr.prevalence +
                               corr.age * age[i] + corr.agesq * age_sq[i] +
                               c(x_cov[i,] %*% corr.beta) +
                               c(z_cov[i,] * corr.alpha) +
                               bias.sim[i,]))
}

stan.data = list(
  n_i = nrow(y.sim),
  n_k = ncol(y.sim),
  x_size = ncol(x_cov),
  x_cov = x_cov,
  z_cov = z_cov,
  age = c(age),
  age_sq = c(age_sq),
  y = y.sim
)



if(model == "Corr_x_z"){
  model.corr.x.z = cmdstan_model("./Stan_Files/Correlated_Stan.stan")
  fit.corr.x.z = model.corr.x.z$sample(data = stan.data,
                                       chains = 3,
                                       parallel_chains = 3,
                                       iter_warmup = 2000,
                                       iter_sampling = 2500,
                                       refresh = 100,
                                       init = 0)
  return = fit.corr.x.z$draws(variables = "prevalence",
                              format = "draws_matrix")
  saveRDS(return, file = paste("./Results/Ukraine_sim_corr_x_z_size_",size.in,"_ind_",ind,".rds", sep = ""))
}else if(model == "Corr_nox_z"){
  model.corr.nox.z = cmdstan_model("./Stan_Files/Correlated_Stan_nox_z.stan")
  fit.corr.nox.z = model.corr.nox.z$sample(data = stan.data,
                                           chains = 3,
                                           parallel_chains = 3,
                                           iter_warmup = 2000,
                                           iter_sampling = 2500,
                                           refresh = 100,
                                           init = 0)
  return = fit.corr.nox.z$draws(variables = "prevalence",
                                format = "draws_matrix")
  saveRDS(return, file = paste("./Results/Ukraine_sim_corr_nox_z_size_",size.in,"_ind_",ind,".rds", sep = ""))
}else if(model == "Corr_x_noz"){
  model.corr.x.noz = cmdstan_model("./Stan_Files/Correlated_Stan_x_noz.stan")
  fit.corr.x.noz = model.corr.x.noz$sample(data = stan.data,
                                           chains = 3,
                                           parallel_chains = 3,
                                           iter_warmup = 2000,
                                           iter_sampling = 2500,
                                           refresh = 100,
                                           init = 0)
  return = fit.corr.x.noz$draws(variables = "prevalence",
                                format = "draws_matrix")
  saveRDS(return, file = paste("./Modeling/Results/Ukraine_sim_corr_x_noz_size_",size.in,"_ind_",ind,".rds", sep = ""))
}else if(model == "Uncorr_x_z"){
  model.uncorr.x.z = cmdstan_model("./Stan_Files/Uncorrelated_Stan_x_z.stan")
  fit.uncorr.x.z = model.uncorr.x.z$sample(data = stan.data,
                                           chains = 3,
                                           parallel_chains = 3,
                                           iter_warmup = 2000,
                                           iter_sampling = 2500,
                                           refresh = 100,
                                           init = 0)
  return = fit.uncorr.x.z$draws(variables = "prevalence",
                                format = "draws_matrix")
  saveRDS(return, file = paste("./Modeling/Results/Ukraine_sim_uncorr_x_z_size_",size.in,"_ind_",ind,".rds", sep = ""))
}else if(model == "Uncorr_nox_noz"){
  model.uncorr.nox.noz = cmdstan_model("./Stan_Files/Uncorrelated_Stan_nox_noz.stan")
  fit.uncorr.nox.noz = model.uncorr.nox.noz$sample(data = stan.data,
                                                   chains = 3,
                                                   parallel_chains = 3,
                                                   iter_warmup = 2000,
                                                   iter_sampling = 2500,
                                                   refresh = 100,
                                                   init = 0)
  return = fit.uncorr.nox.noz$draws(variables = "prevalence",
                                    format = "draws_matrix")
  saveRDS(return, file = paste("./Modeling/Results/Ukraine_sim_uncorr_nox_noz_size_",size.in,"_ind_",ind,".rds", sep = ""))
}else{
  print("Invalid model choice")
}


