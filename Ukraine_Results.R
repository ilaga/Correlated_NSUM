library(ggplot2)
## ggcorrplot and ggpubr called below once


## Load data
y = as.matrix(read.csv("./Data/Ukraine_y.csv", header = F))
x_cov = as.matrix(read.csv("./Data/Ukraine_x_cov.csv", header = F))
x_cov = x_cov[,-c(1, 9)] ## Remove intercept and last column
age = x_cov[,2]
raw.age = age
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
subpop.names = c(names(known), "FSW", "MSW", "MSM", "IDU")

prevalences = known / N
pg1.ind = 1:11
Pg1 = sum(prevalences[pg1.ind])

#############################################################################
#############################################################################
##
## Part 0: Load results
##
#############################################################################
#############################################################################



# Overdispersed Model -----------------------------------------------------
overdisp.betas = matrix(readRDS("./Results/Ukraine_overdispersed_betas.rds"), ncol = 15)
overdisp.alphas = matrix(readRDS("./Results/Ukraine_overdispersed_alphas.rds"), ncol = 8451)
overdisp.omegas = matrix(readRDS("./Results/Ukraine_overdispersed_omegas.rds"), ncol = 15)
n.mc.overdisp = nrow(overdisp.betas)

overdisp.prev = overdisp.betas
for(i in 1:nrow(overdisp.prev)){
  C = log(sum(exp(overdisp.prev[i, pg1.ind])) / Pg1)
  overdisp.prev[i,] = overdisp.prev[i,] - C
}
overdisp.size = exp(overdisp.prev) * N


# Correlated Model --------------------------------------------------------

corr.logdi = readRDS("./Results/Ukraine_correlated_logdi.rds")
corr.prevalence = readRDS("./Results/Ukraine_correlated_prevalence.rds")
corr.bias.in = readRDS("./Results/Ukraine_correlated_bias.rds")
corr.prev.mean.in = readRDS("./Results/Ukraine_correlated_prev_mean.rds")
corr.beta = readRDS("./Results/Ukraine_correlated_beta_par.rds")
corr.alpha = readRDS("./Results/Ukraine_correlated_alpha_par.rds")
corr.age = readRDS("./Results/Ukraine_correlated_age_par.rds")
corr.agesq = readRDS("./Results/Ukraine_correlated_age_sq_par.rds")
corr.tauN = readRDS("./Results/Ukraine_correlated_tau_N.rds")
corr.corr.in = readRDS("./Results/Ukraine_correlated_Corr.rds")
n.mc.corr = nrow(corr.prevalence)

corr.prev = corr.prevalence
for(i in 1:nrow(corr.prev)){
  C = log(sum(exp(corr.prev[i, pg1.ind])) / Pg1)
  corr.prev[i,] = corr.prev[i,] - C
}
corr.size = exp(corr.prev) * N


corr.bias = corr.prev.mean = array(NA, dim = c(n.mc.corr, N.i, N.k))
corr.corr = array(NA, dim = c(n.mc.corr, N.k, N.k))
for(i in 1:n.mc.corr){
  corr.bias[i,,] = matrix(corr.bias.in[i,], nrow = N.i, ncol = N.k)
  corr.prev.mean[i,,] = matrix(corr.prev.mean.in[i,], nrow = N.i, ncol = N.k)
  corr.corr[i,,] = matrix(corr.corr.in[i,], nrow = N.k, ncol = N.k)
  print(i)
}
corr.bias.ave = apply(log(corr.bias), c(2,3), mean)
corr.corr.ave = apply(corr.corr, c(2,3), mean)





# Parameter Estimates -----------------------------------------------------

## Find values for Table 2
colnames(corr.age) = subpop.names
colnames(corr.agesq) = subpop.names
colnames(corr.alpha) = subpop.names

round(colMeans(corr.age), digits = 2)
round(apply(corr.age, 2, quantile, probs = c(0.025, 0.975)), digits = 2)
round(colMeans(corr.agesq), digits = 2)
round(apply(corr.agesq, 2, quantile, probs = c(0.025, 0.975)), digits = 2)
round(colMeans(corr.alpha), digits = 2)
round(apply(corr.alpha, 2, quantile, probs = c(0.025, 0.975)), digits = 2)

## and text in Section 5.2
colnames(corr.beta) = c("Gender", "Nationality", "Employed", "Internet", "Secondary Education", "Vocational Education")
round(colMeans(corr.beta), digits = 2)
round(apply(corr.beta, 2, quantile, probs = c(0.025, 0.975)), digits = 2)

## Create Figure 4
colnames(corr.corr.ave) = rownames(corr.corr.ave) = subpop.names
corr.ordered = ggcorrplot::ggcorrplot(corr.corr.ave, type = "lower",
                       method = "square", hc.order = TRUE,
                       legend.title = "Corr") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 15, face = "bold"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))

corr.unordered = ggcorrplot::ggcorrplot(corr.corr.ave, type = "lower",
                       method = "square", hc.order = FALSE,
                       legend.title = "Corr") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 15, face = "bold"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15))
  
ggsave("./Figures/Ukraine_correlation_ordered.jpg", corr.ordered,
       width = 8, height = 8)
ggsave("./Figures/Ukraine_correlation_unordered.jpg", corr.unordered,
       width = 8, height = 8)
## Only the ordered version is used in the manuscript, but the unordered version
## contains the correlation matrix in the original population ordering


# Get posterior simulations -----------------------------------------------
post.sim.corr = array(NA, dim = c(n.mc.corr, N.i, N.k))
for(ind in 1:n.mc.corr){
  for(k in 1:N.k){
    post.sim.corr[ind,,k] = rpois(N.i, corr.prev.mean[ind,,k] *
                                        corr.bias[ind,,k])
  }
  print(ind)
}

post.sim.overdisp = array(NA, dim = c(n.mc.overdisp, N.i, N.k))
for(ind in 1:n.mc.overdisp){
  for(k in 1:N.k){
    omega.k.m1 = 1/(overdisp.omegas[ind,k] - 1)
    p = omega.k.m1 / (omega.k.m1 + 1)
    xi.i.k = omega.k.m1 * exp(overdisp.alphas[ind,] + overdisp.betas[ind,k])
    post.sim.overdisp[ind,,k] = rnbinom(N.i, xi.i.k, p)
  }
  print(ind)
}



#############################################################################
#############################################################################
##
## Part 1: Subpopulation size estimates
## Figure 2
##
#############################################################################
#############################################################################

corr.upper = corr.lower = rep(NA, 11)
overdisp.upper = overdisp.lower = rep(NA, 11)
for(i in 1:11){
  interval = quantile(corr.size[,i],
                      probs = c(0.025, 0.975), na.rm = T)
  corr.upper[i] = interval[2]
  corr.lower[i] = interval[1]
  
  interval = quantile(overdisp.size[,i], probs = c(0.025, 0.975), na.rm = T)
  overdisp.upper[i] = interval[2]
  overdisp.lower[i] = interval[1]
}


corr.dat = data.frame(est = known - colMeans(corr.size, na.rm = T)[1:11],
                                lower = known - corr.lower,
                                upper = known - corr.upper,
                                model = "Correlated",
                                pop = as.character(1:11),
                                popname = names(known),
                                est.stand = (known - colMeans(corr.size, na.rm = T)[1:11])/known,
                                lower.stand = (known - corr.lower) / known,
                                upper.stand = (known - corr.upper) / known)

overdisp.dat = data.frame(est = known - colMeans(overdisp.size)[1:11],
                       lower = known - overdisp.lower,
                       upper = known - overdisp.upper,
                       model = "Overdispered",
                       pop = as.character(1:11),
                       popname = names(known),
                       est.stand = (known - colMeans(overdisp.size)[1:11])/known,
                       lower.stand = (known - overdisp.lower) / known,
                       upper.stand = (known - overdisp.upper) / known)



dat.all = rbind(corr.dat, overdisp.dat)
dat.all$model = factor(dat.all$model,
                       levels = c("Overdispered", "Correlated"))
subpop.order = as.character(order(known, decreasing = T))

gg.loo = ggplot(data = dat.all) + 
  geom_point(aes(x = popname, y = est.stand, group = model), position=position_dodge(width=0.6), size = 3)+
  geom_errorbar(aes(x = popname, ymin= lower.stand,
                    ymax= upper.stand, lty = model), width=.5, size = 1.4,
                position=position_dodge(width=0.6)) +
  geom_abline(slope = 0, intercept = 0, col = "black", lwd = 1.2) + 
  xlab("Subpopulation") + ylab("Relative Error") +
  scale_linetype_manual(name = "Model", na.translate = F,
                          values = c("dotted", "solid")) +
  scale_x_discrete(limits = names(known)[order(known, decreasing = T)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold", vjust = -10),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 10),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        plot.margin = unit(c(1, 1, 2, 2),"cm"))
gg.loo


## Save Figure 2
ggsave("./Figures/Ukraine_LOO.jpg", gg.loo, width = 15, height = 8)





#############################################################################
#############################################################################
##
## Part 2: Unknown subpopulation sizes
##
#############################################################################
#############################################################################

## Create Table 1
signif(colMeans(corr.size)[12:15], digits = 3)
signif(colMeans(overdisp.size)[12:15], digits = 3)

signif(apply(corr.size[,12:15], 2, quantile, probs = c(0.025, 0.975)), digits = 3)
signif(apply(overdisp.size[,12:15], 2, quantile, probs = c(0.025, 0.975)), digits = 3)









#############################################################################
#############################################################################
##
## Part 3: Prob of responses
## Supplementary Figure 2
##
#############################################################################
#############################################################################


### ppp for prob of responses
resp.vec = list(0, 1, 2, 3, 5, 9, 10, 11, 12:max(y))

resp.dat = matrix(NA, nrow = length(resp.vec), ncol = ncol(y))
for(i in 1:length(resp.vec)){
  resp.dat[i,] = apply(y, 2, function(x){mean(x %in% resp.vec[[i]], na.rm = T)})
}

zero.ppp.corr = array(NA, dim = c(length(resp.vec), n.mc.corr, 15))
zero.ppp.overdisp = array(NA, dim = c(length(resp.vec), n.mc.overdisp, 15))

for(i in 1:length(resp.vec)){
  for(ind in 1:n.mc.corr){
    zero.ppp.corr[i,ind,] = apply(post.sim.corr[ind,,], 2,
                                    function(x){mean(x %in% resp.vec[[i]], na.rm = T)})
  }
  print(i)
}

for(i in 1:length(resp.vec)){
  for(ind in 1:n.mc.overdisp){
    zero.ppp.overdisp[i,ind,] = apply(post.sim.overdisp[ind,,], 2,
                                 function(x){mean(x %in% resp.vec[[i]], na.rm = T)})
  }
  print(i)
}


int.names = c(0, 1, 2, 3, 5, 9, 10, 11, "> 11")

## ggplot
dat.df = matrix(NA, nrow = 0, ncol = 7)
corr.interval = apply(zero.ppp.corr, c(1,3), quantile, probs = c(0.025, 0.975))
overdisp.interval = apply(zero.ppp.overdisp, c(1,3), quantile, probs = c(0.025, 0.975))
names.ind.df = rep(NA, 9)
for(i in 1:9){
  corr.means = colMeans(zero.ppp.corr[i,,])
  overdisp.means = colMeans(zero.ppp.overdisp[i,,])
  if(i == 9){
    dat.df = rbind(dat.df, cbind(resp.dat[i,], corr.means, corr.interval[1,i,],
                                 corr.interval[2,i,], "Correlated",
                                 paste0("Prob(Y ",int.names[i],")"),
                                 c(rep("Known", 11), rep("Unknown", 4))))
    dat.df = rbind(dat.df, cbind(resp.dat[i,], overdisp.means, overdisp.interval[1,i,],
                                 overdisp.interval[2,i,], "Overdispersed",
                                 paste0("Prob(Y ",int.names[i],")"),
                                 c(rep("Known", 11), rep("Unknown", 4))))
    names.ind.df[i] = paste0("P(Y ",int.names[i],")")
  }else{
    dat.df = rbind(dat.df, cbind(resp.dat[i,], corr.means, corr.interval[1,i,],
                                 corr.interval[2,i,], "Correlated",
                                 paste0("Prob(Y = ",int.names[i],")"),
                                 c(rep("Known", 11), rep("Unknown", 4))))
    dat.df = rbind(dat.df, cbind(resp.dat[i,], overdisp.means, overdisp.interval[1,i,],
                                 overdisp.interval[2,i,], "Overdispersed",
                                 paste0("Prob(Y = ",int.names[i],")"),
                                 c(rep("Known", 11), rep("Unknown", 4))))
    names.ind.df[i] = paste0("P(Y = ",int.names[i],")")
  }

}
dat.df = as.data.frame(dat.df)
names(dat.df) = c("Truth", "Val", "Lower", "Upper", "Model", "Ind", "Known")
dat.df$Truth = as.numeric(dat.df$Truth)
dat.df$Val = as.numeric(dat.df$Val)
dat.df$Lower = as.numeric(dat.df$Lower)
dat.df$Upper = as.numeric(dat.df$Upper)
dat.df$Model = as.factor(dat.df$Model)
dat.df$Ind = as.factor(dat.df$Ind)
levels(dat.df$Ind) = names.ind.df
dat.df$Known = as.factor(dat.df$Known)
dat.df$KnownModel = paste0(dat.df$Known, dat.df$Model)
dat.df$KnownModel = factor(dat.df$KnownModel, levels = c("KnownCorrelated",
                                                         "UnknownCorrelated",
                                                         "KnownOverdispersed",
                                                         "UnknownOverdispersed"))

gg.all = ggplot(dat.df, aes(y = Truth, x = Val, shape = KnownModel)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper)) +
  geom_point(size = 3) + 
  facet_wrap(~Ind, scales = "free") +
  geom_abline(slope = 1, intercept = 0, lwd = 1.3) +
  scale_shape_manual(name = "Model", labels = c("Correlated Known",
                                                "Correlated Unknown",
                                                "Overdispersed Known",
                                                "Overdispersed Unknown"),
                       values = c(19, 17, 1, 2)) +
  xlab("Simulated proportion") + ylab("True proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16, face = "bold"))
gg.all

## Save Supplementary Figure 2
ggsave("./Figures/PPP_all_probs.jpg", gg.all,
       width = 20, height = 10)



#############################################################################
#############################################################################
##
## Part 4: ppp Correlation
## Supplementary Figure 3
##
#############################################################################
#############################################################################

dat.cor = cor(y)

ppp.corr.mat = ppp.overdisp.mat = matrix(0, nrow = 15, ncol = 15)
for(i in 1:n.mc.corr){
  cor.z = cor(post.sim.corr[i,,])
  ppp.corr.mat = ppp.corr.mat + (cor.z > dat.cor)
  
  print(i)
}

for(i in 1:n.mc.overdisp){
  cor.overdisp = cor(post.sim.overdisp[i,,])
  ppp.overdisp.mat = ppp.overdisp.mat + (cor.overdisp > dat.cor)
  
  print(i)
}

ppp.corr.mat = ppp.corr.mat/n.mc.corr
ppp.overdisp.mat = ppp.overdisp.mat/n.mc.overdisp

corr.ppp = cbind(ppp.corr.mat[upper.tri(ppp.corr.mat)], "Correlated")
overdisp.ppp = cbind(ppp.overdisp.mat[upper.tri(ppp.overdisp.mat)], "Overdispersed")
all.df = data.frame(rbind(corr.ppp, overdisp.ppp))
names(all.df) = c("Val", "Model")
all.df$Model = factor(all.df$Model, levels = c("Overdispersed", "Correlated"))
all.df$Val = as.numeric(all.df$Val)
all.df$Truth = rep(as.numeric(dat.cor[upper.tri(dat.cor)]), 2)

gg.corr.ppp = ggplot(all.df) + geom_point(aes(x = Val, y = Truth), size = 2) +
  facet_grid(~Model) +
  xlab("Posterior predictive p-value") + ylab("Observed Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 16, face = "bold"))

## Save Supplementary Figure 2
ggsave(paste0("./Figures/Correlation_ppp.jpg"), gg.corr.ppp,
       width = 20, height = 10)

#############################################################################
#############################################################################
##
## Part 5: Surrogate Residuals
## Figure 3
##
#############################################################################
#############################################################################


# Corr Z Model -------------------------------------------------------
## Simulate posterior resid and calculate fay.corr, faym1.corr, etc

# post.sim.corr = array(NA, dim = c(n.mc.corr, N.i, N.k))
post.sim.corr = faym1.corr = fay.corr = faym1.corr.cond =
  fay.corr.cond = array(NA, dim = c(n.mc.corr, N.i, N.k))
for(ind in 1:n.mc.corr){
  for(k in 1:N.k){
    post.sim.corr[ind,,k] = rpois(N.i, corr.prev.mean[ind,,k] *
                                      corr.bias[ind,,k])
    faym1.corr[ind,,k] = ppois(post.sim.corr[ind,,k] - 1, corr.prev.mean[ind,,k] *
                                     corr.bias[ind,,k])
    fay.corr[ind,,k] = ppois(post.sim.corr[ind,,k], corr.prev.mean[ind,,k] *
                                   corr.bias[ind,,k])
    faym1.corr.cond[ind,,k] = ppois(y[,k] - 1, corr.prev.mean[ind,,k] *
                                          corr.bias[ind,,k])
    fay.corr.cond[ind,,k] = ppois(y[,k], corr.prev.mean[ind,,k] *
                                        corr.bias[ind,,k])
  }
  print(ind)
}

## Simulate S.unconditional
S.unconditional.corr = matrix(NA, nrow = N.i, ncol = N.k)
for(i in 1:N.i){
  for(k in 1:N.k){
    S.unconditional.corr[i,k] = mean(runif(n.mc.corr, min = faym1.corr[,i,k], max = fay.corr[,i,k]))
  }
  print(i)
}



## Create plots for Figure 3 (b) and (d)
sample.ind = sample(1:nrow(faym1.corr.cond), 25)
for(k in c(1, 15)){
  S.sim.df = matrix(NA, nrow = 0, ncol = 3)
  for(i in sample.ind){
    S.sim.tmp = runif(N.i, min = faym1.corr.cond[1,,k],
                      max = fay.corr.cond[1,,k])
    S.sim.1 = cbind(S.sim.tmp, S.sim.tmp - S.unconditional.corr[,k], i)
    S.sim.df = rbind(S.sim.df, S.sim.1)
  }
  
  S.sim.df = data.frame(S.sim.df, age = raw.age)
  names(S.sim.df) = c("S.sim", "resid", "sample", "age")
  S.sim.df$sample = factor(S.sim.df$sample)
  
  gg.corr = ggplot(S.sim.df, aes(x = age, y = resid, group = sample)) + geom_point(data = subset(S.sim.df, sample == sample.ind[1])) +
    geom_smooth(lwd = 1.3, se = F) + xlab("Age") + ylab("Surrogate residual") +
    theme(axis.title = element_text(size = 33, face = "bold"),
          axis.text = element_text(size = 28, face = "bold"),
          title = element_text(size = 33, face = "bold")) +
    ggtitle(paste0("Correlated: ",subpop.names[k]))
  
  ggsave(paste0("./Figures/Age_correlated_",subpop.names[k],".jpg"), gg.corr,
         width = 9, height = 9)
}







# Overdispersed Model -------------------------------------------------------
## Simulate posterior resid and calculate fay.overdisp, faym1.overdisp, etc

post.sim.overdisp = faym1.overdisp = fay.overdisp = faym1.overdisp.cond = fay.overdisp.cond =
  array(NA, dim = c(n.mc.overdisp, N.i, N.k))
for(ind in 1:n.mc.overdisp){
  for(k in 1:N.k){
    omega.k.m1 = 1/(overdisp.omegas[ind,k] - 1)
    p = omega.k.m1 / (omega.k.m1 + 1)
    xi.i.k = omega.k.m1 * exp(overdisp.alphas[ind,] + overdisp.betas[ind,k])
    post.sim.overdisp[ind,,k] = rnbinom(N.i, xi.i.k, p)
    faym1.overdisp[ind,,k] = pnbinom(post.sim.overdisp[ind,,k] - 1, xi.i.k, p)
    fay.overdisp[ind,,k] = pnbinom(post.sim.overdisp[ind,,k], xi.i.k, p)
    faym1.overdisp.cond[ind,,k] = pnbinom(y[,k] - 1, xi.i.k, p)
    fay.overdisp.cond[ind,,k] = pnbinom(y[,k], xi.i.k, p)
  }
  print(ind / n.mc.overdisp * 100)
}

## Simulate S.unconditional
S.unconditional.overdisp = matrix(NA, nrow = N.i, ncol = N.k)
for(i in 1:N.i){
  for(k in 1:N.k){
    S.unconditional.overdisp[i,k] = mean(runif(n.mc.overdisp, min = faym1.overdisp[,i,k], max = fay.overdisp[,i,k]))
  }
  print(i)
}


## Create plots for Figure 3 (a) and (c)
sample.ind = sample(1:nrow(faym1.overdisp.cond), 25)
for(k in c(1, 15)){
  S.sim.df = matrix(NA, nrow = 0, ncol = 3)
  for(i in sample.ind){
    S.sim.tmp = runif(N.i, min = faym1.overdisp.cond[1,,k],
                      max = fay.overdisp.cond[1,,k])
    S.sim.1 = cbind(S.sim.tmp, S.sim.tmp - S.unconditional.overdisp[,k], i)
    S.sim.df = rbind(S.sim.df, S.sim.1)
  }
  
  S.sim.df = data.frame(S.sim.df, age = raw.age)
  names(S.sim.df) = c("S.sim", "resid", "sample", "age")
  S.sim.df$sample = factor(S.sim.df$sample)
  
  gg.overdisp = ggplot(S.sim.df, aes(x = age, y = resid, group = sample)) + geom_point(data = subset(S.sim.df, sample == sample.ind[1])) +
    geom_smooth(lwd = 1.6, se = F) + xlab("Age") + ylab("Surrogate residual") +
    theme(axis.title = element_text(size = 33, face = "bold"),
          axis.text = element_text(size = 28, face = "bold"),
          title = element_text(size = 33, face = "bold")) +
    ggtitle(paste0("Overdispersed: ",subpop.names[k]))
  
  ggsave(paste0("./Figures/Age_overdispersed_",subpop.names[k],".jpg"), gg.overdisp,
         width = 9, height = 9)
}







#############################################################################
#############################################################################
##
## Part 5: Conditional mean and std
## Supplementary Figure 1
##
#############################################################################
#############################################################################
post.sim.corr.pos = post.sim.corr
post.sim.corr.pos[post.sim.corr.pos == 0] = NA
post.sim.overdisp.pos = post.sim.overdisp
post.sim.overdisp.pos[post.sim.overdisp.pos == 0] = NA
y.poss = y
y.poss[y.poss == 0] = NA

subpop.order = c(order(known, decreasing = T), 12, 13, 14, 15)
plots = list()
ind = 1
sample.ind = sample(1:n.mc.corr, 500)
for(i in subpop.order){
  plots[[ind]] = local({
    local.i = i
    mean.rep.1 = apply(post.sim.corr.pos[sample.ind,,local.i], 1, mean, na.rm = T)
    sd.rep.1 = apply(post.sim.corr.pos[sample.ind,,local.i], 1, sd, na.rm = T)
    mean.rep.2 = apply(post.sim.overdisp.pos[sample.ind,,local.i], 1, mean, na.rm = T)
    sd.rep.2 = apply(post.sim.overdisp.pos[sample.ind,,local.i], 1, sd, na.rm = T)
    all.df = data.frame(mean = c(mean.rep.1, mean.rep.2),
                        sd = c(sd.rep.1, sd.rep.2),
                        Model = c(rep("Correlated", length(mean.rep.1)),
                                  rep("Overdispersed", length(mean.rep.2))))
    truth.point = data.frame(mean = mean(y.poss[,local.i], na.rm = T),
                             sd = sd(y.poss[,local.i], na.rm = T),
                             Model = "Truth",
                             scale = 3)
    
    corr.df = subset(all.df, Model == "Correlated")
    overdisp.df = subset(all.df, Model == "Overdispersed")

    p = ggplot(all.df) +
      xlab("Mean") + ylab("Std Dev") +
      ggtitle(subpop.names[local.i]) +
      geom_point(data = all.df, mapping = aes(x = mean, y = sd, fill = Model), shape = 21, size = 2.2, stroke = 0.1) +
      geom_point(mapping = aes(x = mean, y = sd, group = Model), data = truth.point, size = 5, col = "black", shape = 18) +
      scale_fill_discrete(name = "Model") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
            axis.text.y = element_text(hjust = 1, size = 16, face = "bold"),
            axis.title.x = element_text(size = 20, face = "bold"),
            axis.title.y = element_text(size = 20, face = "bold"),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 20),
            title = element_text(size = 22, face = "bold"),
            strip.text = element_text(size = 16, face = "bold"))

    print(p)
  })
  
  ind = ind + 1
}


arrange.plot = ggpubr::ggarrange(plots[[1]], plots[[2]], plots[[3]],
                  plots[[4]], plots[[5]], plots[[6]],
                  plots[[7]], plots[[8]], plots[[9]],
                  plots[[10]], plots[[11]], plots[[12]],
                  plots[[13]], plots[[14]], plots[[15]],
                  common.legend = TRUE, legend = "right")

## Save Supplementary Figure 1
ggsave("./Figures/Conditional_stats_scatter.jpg", arrange.plot,
       width = 20, height = 10)
