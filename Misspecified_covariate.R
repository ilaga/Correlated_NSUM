
library(optimx)


## Initialize parameters
k = 3
betas = c(-3, 1, -0.1)
beta.slope = -0.5
## Decompose into tau and Corr for transformation
sigma.sq.list = list(
  rep(0.1, k),
  rep(0.5, k),
  rep(1, k),
  rep(1.5, k)
)
d.sd = sqrt(0.2)

beta.est.mat = beta.est.mat.true = list()
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
  mvt.sd = sqrt(log(1 + tau^2))
  
  
  ###############
  ## Get theoretical probability values
  ###############
  
  p.theta.func.normal.double.truth = function(gamma.in, sigma.in, beta.in, y.in, x.in){
    mvt.mean.tmp = log(1 / sqrt(1 + sigma.in^2))
    mvt.sd.tmp = sqrt(log(1 + sigma.in^2))
    return(dpois(y.in, exp(beta.in + beta.slope * x.in + gamma.in)) *
             dnorm(gamma.in, mean = mvt.mean.tmp, sd = sqrt(mvt.sd.tmp^2 + d.sd^2))) ## Uncorrelated Normal
  }
  
  integrals.vec.x0.1 = integrals.vec.x1.1 = rep(NA, y.max.1 + 1)
  integrals.vec.x0.2 = integrals.vec.x1.2 = rep(NA, y.max.2 + 1)
  integrals.vec.x0.3 = integrals.vec.x1.3 = rep(NA, y.max.3 + 1)
  for(i in 1:(y.max.1 + 1)){
    integrals.vec.x0.1[i] = integrate(function(x) p.theta.func.normal.double.truth(x,
                                                                                   sigma.in = tau[1],
                                                                                   beta.in = betas[1],
                                                                                   y.in = i-1,
                                                                                   x.in = 1),-Inf,Inf)$value
    
    integrals.vec.x1.1[i] = integrate(function(x) p.theta.func.normal.double.truth(x,
                                                                                   sigma.in = tau[1],
                                                                                   beta.in = betas[1],
                                                                                   y.in = i-1,
                                                                                   x.in = -1),-Inf,Inf)$value
  }
  
  for(i in 1:(y.max.2 + 1)){
    integrals.vec.x0.2[i] = integrate(function(x) p.theta.func.normal.double.truth(x,
                                                                                   sigma.in = tau[2],
                                                                                   beta.in = betas[2],
                                                                                   y.in = i-1,
                                                                                   x.in = 1),-Inf,Inf)$value
    
    integrals.vec.x1.2[i] = integrate(function(x) p.theta.func.normal.double.truth(x,
                                                                                   sigma.in = tau[2],
                                                                                   beta.in = betas[2],
                                                                                   y.in = i-1,
                                                                                   x.in = -1),-Inf,Inf)$value
  }
  
  for(i in 1:(y.max.3 + 1)){
    integrals.vec.x0.3[i] = integrate(function(x) p.theta.func.normal.double.truth(x,
                                                                                   sigma.in = tau[3],
                                                                                   beta.in = betas[3],
                                                                                   y.in = i-1,
                                                                                   x.in = 1),-Inf,Inf)$value
    
    integrals.vec.x1.3[i] = integrate(function(x) p.theta.func.normal.double.truth(x,
                                                                                   sigma.in = tau[3],
                                                                                   beta.in = betas[3],
                                                                                   y.in = i-1,
                                                                                   x.in = -1),-Inf,Inf)$value
  }
  
  
  p.theta.pos.x0 = p.theta.pos.x1 = rep(NA, nrow(poss.vals))
  for(i in 1:nrow(poss.vals)){
    p.theta.pos.x0[i] = integrals.vec.x0.1[poss.vals[i,1]+1] *
      integrals.vec.x0.2[poss.vals[i,2]+1] *
      integrals.vec.x0.3[poss.vals[i,3]+1]
    
    p.theta.pos.x1[i] = integrals.vec.x1.1[poss.vals[i,1]+1] *
      integrals.vec.x1.2[poss.vals[i,2]+1] *
      integrals.vec.x1.3[poss.vals[i,3]+1]
  }
  
  
  which.keep = which(p.theta.pos.x0 > 1e-15 & p.theta.pos.x1 > 1e-15)
  p.theta.pos.x0 = p.theta.pos.x0[which.keep]
  p.theta.pos.x1 = p.theta.pos.x1[which.keep]
  y.poss = as.matrix(poss.vals[which.keep,])
  
  
  
  ###############
  ## Step 4: Calculate and maximize expectation
  ###############
  
  p.theta.func.normal.double = function(gamma.in, sigma.in, d.sd.in, beta.slope.in,
                                        beta.in, y.in, x.in){
    mvt.mean.tmp = log(1 / sqrt(1 + sigma.in^2))
    mvt.sd.tmp = sqrt(log(1 + sigma.in^2))
    return(dpois(y.in, exp(beta.in + beta.slope.in * x.in + gamma.in)) *
             dnorm(gamma.in, mean = mvt.mean.tmp, sd = sqrt(mvt.sd.tmp^2 + d.sd.in^2))) ## Uncorrelated Normal
  }
  
  
  
  ##############################################################################
  ##############################################################################
  ## Without X
  
  expectation.func.normal.nox = function(par){
    betas.in = par[1:3]
    sigmas.in = exp(par[4:6])
    d.sd.in = exp(par[7])
    
    for(i in 1:(max(y.poss[,1]) + 1)){
      integrals.vec.x0.1[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[1],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[1],
                                                                               beta.slope.in = 0,
                                                                               y.in = i-1,
                                                                               x.in = 0),-Inf,Inf)$value
    }
    
    
    for(i in 1:(max(y.poss[,2]) + 1)){
      integrals.vec.x0.2[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[2],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[2],
                                                                               beta.slope.in = 0,
                                                                               y.in = i-1,
                                                                               x.in = 0),-Inf,Inf)$value
    }
    
    for(i in 1:(max(y.poss[,3]) + 1)){
      integrals.vec.x0.3[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[3],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[3],
                                                                               beta.slope.in = 0,
                                                                               y.in = i-1,
                                                                               x.in = 0),-Inf,Inf)$value
    }
    
    p.theta.pos.x0.tmp = rep(NA, nrow(y.poss))
    for(i in 1:nrow(y.poss)){
      p.theta.pos.x0.tmp[i] = integrals.vec.x0.1[y.poss[i,1]+1] *
        integrals.vec.x0.2[y.poss[i,2]+1] *
        integrals.vec.x0.3[y.poss[i,3]+1]
    }
    
    
    expect.val = 0
    for(ind in 1:nrow(y.poss)){
      weight.x0 = p.theta.pos.x0[ind]
      weight.x1 = p.theta.pos.x1[ind]
      
      add.val.x0 = log(p.theta.pos.x0.tmp[ind])
      expect.val = expect.val + add.val.x0 * weight.x0 + add.val.x0 * weight.x1
    }
    
    print(expect.val)
    print(par)
    
    return(expect.val)
  }
  
  
  ##############################################################################
  ##############################################################################
  ## With X
  expectation.func.normal.x = function(par){
    betas.in = par[1:3]
    sigmas.in = exp(par[4:6])
    beta.slope.in = par[7]
    d.sd.in = exp(par[8])
    
    integrals.vec.x0.1 = integrals.vec.x0.2 = integrals.vec.x0.3 = rep(NA, max(y.poss[,1]) + 1)
    integrals.vec.x1.1 = integrals.vec.x1.2 = integrals.vec.x1.3 = rep(NA, max(y.poss[,2]) + 1)
    for(i in 1:(max(y.poss[,1]) + 1)){
      integrals.vec.x0.1[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[1],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[1],
                                                                               beta.slope.in = beta.slope.in,
                                                                               y.in = i-1,
                                                                               x.in = 1),-Inf,Inf)$value
    }
    
    for(i in 1:(max(y.poss[,2]) + 1)){
      integrals.vec.x0.2[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[2],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[2],
                                                                               beta.slope.in = beta.slope.in,
                                                                               y.in = i-1,
                                                                               x.in = 1),-Inf,Inf)$value
    }
    
    for(i in 1:(max(y.poss[,3]) + 1)){
      integrals.vec.x0.3[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[3],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[3],
                                                                               beta.slope.in = beta.slope.in,
                                                                               y.in = i-1,
                                                                               x.in = 1),-Inf,Inf)$value
    }
    
    
    for(i in 1:(max(y.poss[,1]) + 1)){
      integrals.vec.x1.1[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[1],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[1],
                                                                               beta.slope.in = beta.slope.in,
                                                                               y.in = i-1,
                                                                               x.in = -1),-Inf,Inf)$value
    }
    
    for(i in 1:(max(y.poss[,2]) + 1)){
      integrals.vec.x1.2[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[2],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[2],
                                                                               beta.slope.in = beta.slope.in,
                                                                               y.in = i-1,
                                                                               x.in = -1),-Inf,Inf)$value
    }
    
    for(i in 1:(max(y.poss[,3]) + 1)){
      integrals.vec.x1.3[i] = integrate(function(x) p.theta.func.normal.double(x,
                                                                               sigma.in = sigmas.in[3],
                                                                               d.sd.in = d.sd.in,
                                                                               beta.in = betas.in[3],
                                                                               beta.slope.in = beta.slope.in,
                                                                               y.in = i-1,
                                                                               x.in = -1),-Inf,Inf)$value
    }
    
    p.theta.pos.x0.tmp = p.theta.pos.x1.tmp = rep(NA, nrow(y.poss))
    for(i in 1:nrow(y.poss)){
      p.theta.pos.x0.tmp[i] = integrals.vec.x0.1[y.poss[i,1]+1] *
        integrals.vec.x0.2[y.poss[i,2]+1] *
        integrals.vec.x0.3[y.poss[i,3]+1]
      
      p.theta.pos.x1.tmp[i] = integrals.vec.x1.1[y.poss[i,1]+1] *
        integrals.vec.x1.2[y.poss[i,2]+1] *
        integrals.vec.x1.3[y.poss[i,3]+1]
    }
    
    
    expect.val = 0
    for(ind in 1:nrow(y.poss)){
      weight.x0 = p.theta.pos.x0[ind]
      weight.x1 = p.theta.pos.x1[ind]
      
      add.val.x0 = log(p.theta.pos.x0.tmp[ind])
      add.val.x1 = log(p.theta.pos.x1.tmp[ind])
      expect.val = expect.val + add.val.x0 * weight.x0 + add.val.x1 * weight.x1
    }
    
    print(expect.val)
    print(par)
    
    return(expect.val)
  }
  
  optim.normal.nox = optimr(par = c(betas, log(tau), log(d.sd)), fn = expectation.func.normal.nox,
                            control = list(fnscale = -1), method = "Nelder-Mead")
  
  optim.normal.x = optimr(par = c(betas, log(tau), beta.slope, log(d.sd)), fn = expectation.func.normal.x,
                            control = list(fnscale = -1), method = "Nelder-Mead")
  
  beta.est.mat[[loop.ind]] = c(optim.normal.nox$par[1:3], exp(optim.normal.nox$par[4:6]))
  beta.est.mat.true[[loop.ind]] = c(optim.normal.x$par[1:3], exp(optim.normal.x$par[4:6]))
  
  
  print(c("Done with overall loop:", loop.ind))
}


## Values for Supplementary Table 2

round(100 * (c(betas, sqrt(sigma.sq.list[[1]])) - beta.est.mat.true[[1]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[2]])) - beta.est.mat.true[[2]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[3]])) - beta.est.mat.true[[3]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[4]])) - beta.est.mat.true[[4]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)

round(100 * (c(betas, sqrt(sigma.sq.list[[1]])) - beta.est.mat[[1]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[2]])) - beta.est.mat[[2]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[3]])) - beta.est.mat[[3]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)
round(100 * (c(betas, sqrt(sigma.sq.list[[4]])) - beta.est.mat[[4]]) / c(betas, sqrt(sigma.sq.list[[1]])), digits = 0)


