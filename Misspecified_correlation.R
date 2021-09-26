
library(mvtnorm)
library(optimx)
library(cubature)


## Initialize parameters
k = 3
betas = c(-3, 1, -0.1)
Corr = rbind(c(1, 2/3, 0),
             c(2/3, 1, 0),
             c(0, 0, 1))
## Decompose into tau and Corr for transformation
sigma.sq.list = list(
  rep(0.1, k),
  rep(0.5, k),
  rep(1, k),
  rep(1.5, k)
)
d.sd = sqrt(0.2)

beta.est.mat = list()
variance.mat = list()
model.variance.mat = list()

y.max.1 = 32
y.max.2 = 200
y.max.3 = 80
poss.vals = as.data.frame(expand.grid(0:y.max.1, 0:y.max.2, 0:y.max.3))
for(loop.ind in 1:length(sigma.sq.list)){
  beta.est.tmp = matrix(NA, nrow = 3, ncol = 2)
  sigma.sq = sigma.sq.list[[loop.ind]]
  
  tau = sqrt(sigma.sq)
  mvt.mean = log(1 / sqrt(1 + tau^2))
  mvt.tau = sqrt(log(1 + tau^2))
  Sigma.mat = diag(mvt.tau) %*% Corr %*% diag(mvt.tau)
  
  ###############
  ## Get theoretical probability values
  ###############
  
  p.theta.func.normal = function(gamma.in, y.in){
    return(prod(dpois(y.in, exp(betas + gamma.in))) *
             dmvnorm(gamma.in, mean = mvt.mean, sigma = Sigma.mat + diag(rep(d.sd^2, 3)))) ## Correlated Normal
  }
  
  start = Sys.time()
  p.theta.pos.int = rep(NA, nrow(poss.vals))
  for(i in 1:nrow(poss.vals)){
    integrals.eval.1 = cubintegrate(p.theta.func.normal, lower = rep(-Inf, 3),
                                    upper = rep(Inf, 3),
                                    y.in = as.numeric(poss.vals[i,]),
                                    method = "hcubature")$integral
    
    p.theta.pos.int[i] = integrals.eval.1
    print(i / length(p.theta.pos.int) * 100)
  }
  Sys.time() - start
  
  
  which.keep = which(p.theta.pos.int > 1e-10)
  p.theta.pos = p.theta.pos.int[which.keep]
  y.poss = as.matrix(poss.vals[which.keep,])
  
  
  
  
  
  
  ###############
  ## Step 4: Calculate and maximize expectation
  ###############
  
  p.theta.func.corr.fixed = function(gamma.in, sigma.in, d.sd.in, betas.in, y.in){
    mvt.mean.tmp = log(1 / sqrt(1 + sigma.in^2))
    mvt.sd.tmp = sqrt(log(1 + sigma.in^2))
    return(prod(dpois(y.in, exp(betas.in + gamma.in))) *
             dmvnorm(gamma.in, mean = mvt.mean.tmp, sigma = diag(mvt.sd.tmp) %*% Corr %*% diag(mvt.sd.tmp) + diag(rep(d.sd.in^2, 3)))) ## Correlated Normal
  }
  
  
  p.theta.func.normal = function(gamma.in, sigma.in,
                                 beta.in, y.in, d.sd.in){
    mvt.mean.tmp = log(1 / sqrt(1 + sigma.in^2))
    mvt.sd.tmp = sqrt(log(1 + sigma.in^2))
    return(dpois(y.in, exp(beta.in + gamma.in)) *
             dnorm(gamma.in, mean = mvt.mean.tmp, sd = sqrt(mvt.sd.tmp^2 + d.sd.in^2))) ## Uncorrelated Normal
  }
  
  ##############################################################################
  ##############################################################################
  ## Normal
  integrals.eval.1 = rep(NA, nrow(y.poss))
  
  expectation.corr.fixed = function(par){
    betas.in = par[1:3]
    sigma.in = exp(par[4:6])
    d.sd.in = exp(par[7])
    
    p.theta.pos.x.tmp = rep(NA, nrow(y.poss))
    for(i in 1:nrow(y.poss)){
      integrals.eval.1 = cubintegrate(p.theta.func.corr.fixed, lower = rep(-Inf, 3),
                                      upper = rep(Inf, 3),
                                      sigma.in = sigma.in,
                                      y.in = as.numeric(y.poss[i,]),
                                      d.sd.in = d.sd.in,
                                      betas.in = betas.in,
                                      method = "hcubature")$integral
      p.theta.pos.x.tmp[i] = log(integrals.eval.1)
    }
    
    expect.val = 0
    for(ind in 1:nrow(y.poss)){
      weight.x0 = p.theta.pos[ind]
      
      add.val.x0 = p.theta.pos.x.tmp[ind]
      expect.val = expect.val + add.val.x0 * weight.x0
    }
    
    print(expect.val)
    print(par)
    
    return(expect.val)
  }
  
  
  expectation.func.normal = function(par){
    betas.in = par[1:3]
    sigmas.in = exp(par[4:6])
    d.sd.in = exp(par[7])
    
    int.tmp.1 = rep(NA, y.max.1 + 1)
    int.tmp.2 = rep(NA, y.max.2 + 1)
    int.tmp.3 = rep(NA, y.max.3 + 1)
    for(i in 1:(y.max.1 + 1)){
      int.tmp.1[i] = integrate(function(x) p.theta.func.normal(x,
                                                               sigma.in = sigmas.in[1],
                                                               d.sd.in = d.sd.in,
                                                               beta.in = betas.in[1],
                                                               y.in = i-1),-Inf,Inf)$value
      
    }
    
    for(i in 1:(y.max.2 + 1)){
      int.tmp.2[i] = integrate(function(x) p.theta.func.normal(x,
                                                               sigma.in = sigmas.in[2],
                                                               d.sd.in = d.sd.in,
                                                               beta.in = betas.in[2],
                                                               y.in = i-1),-Inf,Inf)$value
    }
    
    for(i in 1:(y.max.3 + 1)){
      int.tmp.3[i] = integrate(function(x) p.theta.func.normal(x,
                                                               sigma.in = sigmas.in[3],
                                                               d.sd.in = d.sd.in,
                                                               beta.in = betas.in[3],
                                                               y.in = i-1),-Inf,Inf)$value
    }
    
    p.theta.pos.x.tmp = rep(NA, nrow(y.poss))
    for(i in 1:nrow(y.poss)){
      p.theta.pos.x.tmp[i] = log(int.tmp.1[y.poss[i,1]+1]) +
        log(int.tmp.2[y.poss[i,2]+1]) +
        log(int.tmp.3[y.poss[i,3]+1])
    }
    
    expect.val = 0
    for(ind in 1:nrow(y.poss)){
      weight.x0 = p.theta.pos[ind]
      
      add.val.x0 = p.theta.pos.x.tmp[ind]
      expect.val = expect.val + add.val.x0 * weight.x0
    }
    
    print(expect.val)
    print(par)
    
    return(expect.val)
  }
  
  optim.normal = optimr(par = c(betas, log(tau), log(d.sd)) + rnorm(7, sd = 0.1), fn = expectation.func.normal,
                        control = list(fnscale = -1, reltol = 1e-14), method = "Nelder-Mead")
  
  optim.corr.fixed = optimr(par = c(betas, log(tau), log(d.sd)) + rnorm(7, sd = 0.1), fn = expectation.corr.fixed,
                            control = list(fnscale = -1, reltol = 1e-14), method = "Nelder-Mead")
  
  
  beta.est.mat[[loop.ind]] = c(optim.normal$par[1:3], exp(optim.normal$par[4:6]))
  beta.est.mat.true[[loop.ind]] = c(optim.corr.fixed$par[1:3], exp(optim.corr.fixed$par[4:6]))
  

  print(c("Done with overall loop:", loop.ind))
}


## Values for Supplementary Table 3

round(100 * (c(betas, sqrt(sigma.sq.list[[1]])) - beta.est.mat.true[[1]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[2]])) - beta.est.mat.true[[2]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[3]])) - beta.est.mat.true[[3]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[4]])) - beta.est.mat.true[[4]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)

round(100 * (c(betas, sqrt(sigma.sq.list[[1]])) - beta.est.mat[[1]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[2]])) - beta.est.mat[[2]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[3]])) - beta.est.mat[[3]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[4]])) - beta.est.mat[[4]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)