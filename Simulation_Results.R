
library(ggplot2)
## Load data
y = as.matrix(read.csv("./Data/Ukraine_y.csv", header = F))
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
known = known[keep.known]
n.known = length(known)
##
##########
N.i = nrow(y)
N.k = ncol(y)
known.orig = known

prevalences.sim = colMeans(readRDS("./Results/Ukraine_correlated_prevalence.rds"))
known = prevalences.sim[1:15]
names(known) = c(names(known.orig), "FSW", "MSW", "MSM", "IDU")

dat.all.sizes = matrix(NA, nrow = 0, ncol = 10)
cov.list = list()


size.ind = 1
for(size.in in c(100,300,1000,3000,10000)){
  uncorr.nox.noz.lower = uncorr.nox.noz.upper = uncorr.nox.noz.mean = uncorr.nox.noz.cov = matrix(NA, nrow = 100, ncol = 15)
  uncorr.x.z.lower = uncorr.x.z.upper = uncorr.x.z.mean = uncorr.x.z.cov = matrix(NA, nrow = 100, ncol = 15)
  corr.x.noz.lower = corr.x.noz.upper = corr.x.noz.mean = corr.x.noz.cov = matrix(NA, nrow = 100, ncol = 15)
  corr.nox.z.lower = corr.nox.z.upper = corr.nox.z.mean = corr.nox.z.cov = matrix(NA, nrow = 100, ncol = 15)
  corr.x.z.lower = corr.x.z.upper = corr.x.z.mean = corr.x.z.cov = matrix(NA, nrow = 100, ncol = 15)
  
  uncorr.nox.noz.size.mean = uncorr.x.z.size.mean = corr.x.noz.size.mean = corr.nox.z.size.mean =
    corr.x.z.size.mean = matrix(NA, nrow = 100, ncol = 15)
  uncorr.nox.noz.width = uncorr.x.z.width = corr.x.noz.width = corr.nox.z.width =
    corr.x.z.width = matrix(NA, nrow = 100, ncol = 15)

  for(ind in 1:100){
    uncorr.nox.noz.res = try(readRDS(paste0("./Results/Ukraine_sim_uncorr_nox_noz_size_", size.in, "_ind_",ind,".rds")))
    uncorr.x.z.res = try(readRDS(paste0("./Results/Ukraine_sim_uncorr_x_z_size_", size.in, "_ind_",ind,".rds")))
    corr.nox.z.res = try(readRDS(paste0("./Results/Ukraine_sim_corr_nox_z_size_", size.in, "_ind_",ind,".rds")))
    corr.x.noz.res = try(readRDS(paste0("./Results/Ukraine_sim_corr_x_noz_size_", size.in, "_ind_",ind,".rds")))
    corr.x.z.res = try(readRDS(paste0("./Results/Ukraine_sim_corr_x_z_size_", size.in, "_ind_",ind,".rds")))
    
    ## The following lines are useful if you don't have all 100 simulations
    if("try-error" %in% class(uncorr.nox.noz.res)){
      uncorr.nox.noz.res = matrix(NA, nrow = 7500, ncol = 15)
    }
    if("try-error" %in% class(uncorr.x.z.res)){
      uncorr.x.z.res = matrix(NA, nrow = 7500, ncol = 15)
    }
    if("try-error" %in% class(corr.nox.z.res)){
      corr.nox.z.res = matrix(NA, nrow = 7500, ncol = 15)
    }
    if("try-error" %in% class(corr.x.noz.res)){
      corr.x.noz.res = matrix(NA, nrow = 7500, ncol = 15)
    }
    if("try-error" %in% class(corr.x.z.res)){
      corr.x.z.res = matrix(NA, nrow = 7500, ncol = 15)
    }

    uncorr.nox.noz.mean[ind,] = colMeans(uncorr.nox.noz.res)
    uncorr.x.z.mean[ind,] = colMeans(uncorr.x.z.res)
    corr.nox.z.mean[ind,] = colMeans(corr.nox.z.res)
    corr.x.noz.mean[ind,] = colMeans(corr.x.noz.res)
    corr.x.z.mean[ind,] = colMeans(corr.x.z.res)
    
    for(i in 1:15){
      interval = quantile(uncorr.nox.noz.res[,i],
                          probs = c(0.025, 0.975), na.rm = T)
      uncorr.nox.noz.upper[ind,i] = interval[2]
      uncorr.nox.noz.lower[ind,i] = interval[1]
      uncorr.nox.noz.cov[ind,i] = (known[i] > uncorr.nox.noz.lower[ind,i] &
                                     known[i] < uncorr.nox.noz.upper[ind,i])
      uncorr.nox.noz.width[ind,i] = interval[2] - interval[1]
      
      
      interval = quantile(uncorr.x.z.res[,i],
                          probs = c(0.025, 0.975), na.rm = T)
      uncorr.x.z.upper[ind,i] = interval[2]
      uncorr.x.z.lower[ind,i] = interval[1]
      uncorr.x.z.cov[ind,i] = (known[i] > uncorr.x.z.lower[ind,i] &
                                 known[i] < uncorr.x.z.upper[ind,i])
      uncorr.x.z.width[ind,i] = interval[2] - interval[1]
      
      interval = quantile(corr.nox.z.res[,i],
                          probs = c(0.025, 0.975), na.rm = T)
      corr.nox.z.upper[ind,i] = interval[2]
      corr.nox.z.lower[ind,i] = interval[1]
      corr.nox.z.cov[ind,i] = (known[i] > corr.nox.z.lower[ind,i] &
                                 known[i] < corr.nox.z.upper[ind,i])
      corr.nox.z.width[ind,i] = interval[2] - interval[1]
      
      interval = quantile(corr.x.noz.res[,i],
                          probs = c(0.025, 0.975), na.rm = T)
      corr.x.noz.upper[ind,i] = interval[2]
      corr.x.noz.lower[ind,i] = interval[1]
      corr.x.noz.cov[ind,i] = (known[i] > corr.x.noz.lower[ind,i] &
                                 known[i] < corr.x.noz.upper[ind,i])
      corr.x.noz.width[ind,i] = interval[2] - interval[1]
      
      interval = quantile(corr.x.z.res[,i],
                          probs = c(0.025, 0.975), na.rm = T)
      corr.x.z.upper[ind,i] = interval[2]
      corr.x.z.lower[ind,i] = interval[1]
      corr.x.z.cov[ind,i] = (known[i] > corr.x.z.lower[ind,i] &
                               known[i] < corr.x.z.upper[ind,i])
      corr.x.z.width[ind,i] = interval[2] - interval[1]
    }
  }
  
  
  uncorr.nox.noz.dat = data.frame(model = "Uncorrelated No x, No Z",
                                  X = "No",
                                  Z = "No",
                                  Corr = "No",
                                  Dist = "Normal",
                                  pop = as.character(1:15),
                                  popname = names(known),
                                  est.stand = (colMeans(uncorr.nox.noz.mean, na.rm = T)[1:15]),
                                  lower.stand = apply(uncorr.nox.noz.mean, 2, quantile, probs = c(0.025), na.rm = T),
                                  upper.stand = apply(uncorr.nox.noz.mean, 2, quantile, probs = c(0.975), na.rm = T),
                                  width = colMeans(uncorr.nox.noz.width, na.rm = T),
                                  lower.credible = colMeans(uncorr.nox.noz.lower, na.rm = T),
                                  upper.credible = colMeans(uncorr.nox.noz.upper, na.rm = T),
                                  truth = known)
  
  uncorr.x.z.dat = data.frame(model = "Uncorrelated Yes x, Yes Z",
                              X = "Yes",
                              Z = "Yes",
                              Corr = "No",
                              Dist = "Normal",
                              pop = as.character(1:15),
                              popname = names(known),
                              est.stand = (colMeans(uncorr.x.z.mean, na.rm = T)[1:15]),
                              lower.stand = apply(uncorr.x.z.mean, 2, quantile, probs = c(0.025), na.rm = T),
                              upper.stand = apply(uncorr.x.z.mean, 2, quantile, probs = c(0.975), na.rm = T),
                              width = colMeans(uncorr.x.z.width, na.rm = T),
                              lower.credible = colMeans(uncorr.x.z.lower, na.rm = T),
                              upper.credible = colMeans(uncorr.x.z.upper, na.rm = T),
                              truth = known)
  
  corr.x.noz.dat = data.frame(model = "Correlated Yes x, No Z",
                              X = "Yes",
                              Z = "No",
                              Corr = "Yes",
                              Dist = "Normal",
                              pop = as.character(1:15),
                              popname = names(known),
                              est.stand = (colMeans(corr.x.noz.mean, na.rm = T)[1:15]),
                              lower.stand = apply(corr.x.noz.mean, 2, quantile, probs = c(0.025), na.rm = T),
                              upper.stand = apply(corr.x.noz.mean, 2, quantile, probs = c(0.975), na.rm = T),
                              width = colMeans(corr.x.noz.width, na.rm = T),
                              lower.credible = colMeans(corr.x.noz.lower, na.rm = T),
                              upper.credible = colMeans(corr.x.noz.upper, na.rm = T),
                              truth = known)
  
  corr.nox.z.dat = data.frame(model = "Correlated No x, Yes Z",
                              X = "No",
                              Z = "Yes",
                              Corr = "Yes",
                              Dist = "Normal",
                              pop = as.character(1:15),
                              popname = names(known),
                              est.stand = (colMeans(corr.nox.z.mean, na.rm = T)[1:15]),
                              lower.stand = apply(corr.nox.z.mean, 2, quantile, probs = c(0.025), na.rm = T),
                              upper.stand = apply(corr.nox.z.mean, 2, quantile, probs = c(0.975), na.rm = T),
                              width = colMeans(corr.nox.z.width, na.rm = T),
                              lower.credible = colMeans(corr.nox.z.lower, na.rm = T),
                              upper.credible = colMeans(corr.nox.z.upper, na.rm = T),
                              truth = known)
  
  corr.x.z.dat = data.frame(model = "Correlated Yes x, Yes Z",
                            X = "Yes",
                            Z = "Yes",
                            Corr = "Yes",
                            Dist = "Normal",
                            pop = as.character(1:15),
                            popname = names(known),
                            est.stand = (colMeans(corr.x.z.mean, na.rm = T)[1:15]),
                            lower.stand = apply(corr.x.z.mean, 2, quantile, probs = c(0.025), na.rm = T),
                            upper.stand = apply(corr.x.z.mean, 2, quantile, probs = c(0.975), na.rm = T),
                            width = colMeans(corr.x.z.width, na.rm = T),
                            lower.credible = colMeans(corr.x.z.lower, na.rm = T),
                            upper.credible = colMeans(corr.x.z.upper, na.rm = T),
                            truth = known)

  dat.all = rbind(uncorr.nox.noz.dat,
                  uncorr.x.z.dat,
                  corr.x.noz.dat,
                  corr.nox.z.dat,
                  corr.x.z.dat)
  
  dat.all$model = factor(dat.all$model,
                         levels = c("Uncorrelated No x, No Z",
                                    "Uncorrelated Yes x, Yes Z",
                                    "Correlated Yes x, No Z",
                                    "Correlated No x, Yes Z",
                                    "Correlated Yes x, Yes Z"))
  dat.all$size = size.in
  dat.all.sizes = rbind(dat.all.sizes, dat.all)
  cov.list[[size.ind]] = rbind(colMeans(uncorr.nox.noz.cov, na.rm = T),
                        colMeans(uncorr.x.z.cov, na.rm = T),
                        colMeans(corr.x.noz.cov, na.rm = T),
                        colMeans(corr.nox.z.cov, na.rm = T),
                        colMeans(corr.x.z.cov, na.rm = T))
  size.ind = size.ind + 1
}

## Create Figure 5
dat.all.sizes$size = factor(dat.all.sizes$size, levels = c("100", "300", "1000",
                                               "3000", "10000"))
dat.all.sizes$Dist = factor(dat.all.sizes$Dist, levels = c("Gamma", "Normal"))
dat.all.sizes$popname = factor(dat.all.sizes$popname, levels = names(known))
## Only plot 9 of the subpopulations
dat.subset.sizes = subset(dat.all.sizes, pop %in% c(1, 4, 5, 7, 9, 12, 13, 14, 15))

beta.gg = ggplot(data = dat.subset.sizes) + 
  geom_point(aes(x = size, y = est.stand, col = model, group = model),
             position=position_dodge(width=0.6), size = 3)+
  geom_errorbar(aes(x = size, ymin= lower.stand,
                    ymax= upper.stand, col = model, group = model),
                width=.4, size = 0.8,
                position=position_dodge(width=0.6)) +
  xlab("Sample Size") + ylab("Estimate") +
  scale_color_discrete(labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 16, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 18, face = "bold")) +
  guides(color = guide_legend(override.aes = list(shape = 21),
                              title = "Model")) +
  facet_wrap(~popname, scales = "free_y", nrow = 3) +
  geom_hline(aes(yintercept = truth), lwd = 1.2)
beta.gg
ggsave("./Figures/Simulation_rho.png", beta.gg, dpi = 300, width = 20, height = 10)