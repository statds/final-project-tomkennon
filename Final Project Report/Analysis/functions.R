###############
###Functions###
###############



### gevr_fit_nloc function
## This function fits the gevR model to the data.
## Implemented new nonlinear location parameter capability
source("gevr_fit_nloc.R")




#Gompertz function evaluates parameters using gompertz distribution
gompertz <- function(time, theta){
  theta[1] + theta[2]*exp(-theta[3]*exp(-theta[4]*time))
}


### Log likelihood function for gevr distribution using gompertz curve as location parameter
gevrloglik <- function(theta, data){
  time <- 1:NROW(data)
  mu <- gompertz(time, theta[1:4])
  #if (theta[1] + theta[2] - theta[5]/theta[6] > 1) return(-Inf)
  if (theta[6] < 0 & (theta[1] + theta[2] - theta[5]/theta[6] > 1)) return(-Inf)
  if (theta[6] > 0 & (theta[1] + theta[2] - theta[5]/theta[6] < 0)) return(-Inf)
  if (theta[5] <= 0) return(-Inf)
  log.density <- sum(dgevr(data, loc = mu, scale = theta[5], shape = theta[6], log.d=TRUE))
  return(log.density)
}

gevrloglik_const <- function(theta, data){
  mu <- mean(data[1])
  #if (theta[3] < 0 & (theta[1] - theta[2]/theta[3] > 1)) return(-Inf)
  #if (theta[3] > 0 & (theta[1] - theta[2]/theta[3] < 0)) return(-Inf)
  if (theta[2] <= 0) return(-Inf)
  log.density <- sum(dgevr(data, loc = mu, scale = theta[2], shape = theta[3], log.d=TRUE))
  return(log.density)
}





## Simulates data a given number of bootstraps for a given r order value
## Then evaluates the log likelihood of maximum likelihood evaluated parameters
## The resulting log likelihoods are stored in a vector called "stats"
sim_loglik <- function(r, boots, data, obs_par) {
  stats <- rep(NA,boots)
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=length(data), r = r, loc=gompertz(1:length(data),obs_par[1:4]), scale=obs_par[5], shape=obs_par[6])
    mle_obj <- matrix(nrow=length(data),ncol=r)
    mle_obj <- optim(obs_par,gevrloglik, data=sim_data, control=list(fnscale=-1, trace=FALSE, maxit=999))
    mle_par <- mle_obj$par
    test <- gevrloglik(mle_par, data=sim_data)
    stats[i] <- test
  }
  return(stats)
}



sim_param <- function(r, boots, data, obs_par) {
  stats <- matrix(nrow=boots, ncol=(length(obs_par)+3))
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=length(data), r = r, loc=gompertz(1:length(data),obs_par[1:4]),
                      scale=obs_par[5], shape=obs_par[6])
    mle_obj <- matrix(nrow=length(data),ncol=r)
    mle_obj <- optim(obs_par,gevrloglik, data=sim_data,
                       control=list(fnscale=-1, trace=FALSE, maxit=999))
    conv <- mle_obj$convergence
    mle_par <- mle_obj$par
    mle_val <- mle_obj$value
    time = 1:length(data)
    mu <- gompertz(time, mle_par[1:4])
    non_loc_data = data - mu
    ks = ks.test(non_loc_data,pgev,scale=mle_par[5],shape=mle_par[6])
    stats[i,] <- c(mle_par,mle_val,ks$statistic,conv)
  }
  return(stats)
}

sim_param2 <- function(r, boots, data, obs_par) {
  stats <- matrix(nrow=boots, ncol=(length(obs_par)+3))
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=nrow(data), r = r, loc=gompertz(1:nrow(data),obs_par[1:4]),
                      scale=obs_par[5], shape=obs_par[6])
    mle_obj <- matrix(nrow=nrow(data),ncol=r)
    mle_obj <- optim(obs_par,gevrloglik, data=sim_data,
                     control=list(fnscale=-1, trace=FALSE, maxit=999))
    conv <- mle_obj$convergence
    mle_par <- mle_obj$par
    mle_val <- mle_obj$value
    time = 1:nrow(data)
    mu <- gompertz(time, mle_par[1:4])
    non_loc_data = data - mu
    ks = ks.test(non_loc_data,pgev,scale=mle_par[5],shape=mle_par[6])
    stats[i,] <- c(mle_par,mle_val,ks$statistic,conv)
  }
  return(stats)
}

sim_param_constant <- function(r, boots, data, obs_par) {
  stats <- matrix(nrow=boots, ncol=(length(obs_par)+2))
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=length(data), r = r, loc=obs_par[1],
                      scale=obs_par[2], shape=obs_par[3])
    mle_obj <- gevrFit(sim_data, method = "mle")
    mle_par <- mle_obj$par.ests
    mle_val <- mle_obj$nllh.final
    mu <- mle_par[1]
    non_loc_data = data - mu
    ks = ks.test(non_loc_data,pgev,scale=mle_par[2],shape=mle_par[3])
    stats[i,] <- c(mle_par,mle_val,ks$statistic)
  }
  return(stats)
}

sim_param_constant2 <- function(r, boots, data, obs_par) {
  stats <- matrix(nrow=boots, ncol=(length(obs_par)+3))
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=length(data), r = r, loc=obs_par[1],
                      scale=obs_par[2], shape=obs_par[3])
    mle_obj <- optim(obs_par, gevrloglik_const, data=sim_data,
                     control=list(fnscale=-1, trace=TRUE, maxit=999))
    conv <- mle_obj$convergence
    mle_par <- mle_obj$par
    mle_val <- mle_obj$value
    mu <- mle_par[1]
    non_loc_data = data - mu
    ks = ks.test(non_loc_data,pgev,scale=mle_par[2],shape=mle_par[3])
    stats[i,] <- c(mle_par,mle_val,ks$statistic,conv)
  }
  return(stats)
}


sim_param_constant3 <- function(r, boots, data, obs_par) {
  stats <- matrix(nrow=boots, ncol=(length(obs_par)+3))
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=nrow(data), r = r, loc=obs_par[1],
                      scale=obs_par[2], shape=obs_par[3])
    mle_obj <- optim(obs_par, gevrloglik_const, data=sim_data,
                     control=list(fnscale=-1, trace=TRUE, maxit=999))
    conv <- mle_obj$convergence
    mle_par <- mle_obj$par
    mle_val <- mle_obj$value
    mu <- mle_par[1]
    non_loc_data = data - mu
    ks = ks.test(non_loc_data,pgev,scale=mle_par[2],shape=mle_par[3])
    stats[i,] <- c(mle_par,mle_val,ks$statistic,conv)
  }
  return(stats)
}


sim_param_constant4 <- function(r, boots, data, obs_par) {
  stats <- matrix(nrow=boots, ncol=(length(obs_par)+2))
  for(i in 1:boots){
    sim_data <- rep(NA, boots)
    sim_data <- rgevr(n=length(data), r = r, loc=obs_par[1],
                      scale=obs_par[2], shape=obs_par[3])
    mle_obj <- gevrFit(sim_data, method = "mle")
    mle_par <- mle_obj$par.ests
    mle_val <- mle_obj$nllh.final
    mu <- mle_par[1]
    non_loc_data = data - mu
    ks = ks.test(non_loc_data,pgev,scale=mle_par[2],shape=mle_par[3])
    stats[i,] <- c(mle_par,mle_val,ks$statistic)
  }
  return(stats)
}
