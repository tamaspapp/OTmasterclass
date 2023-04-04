source("samplers.R")
source("helpers.R")

#########
# Setup #
#########

# Target log-density
logpi <- function(x){-0.5*sum(x^2)}

# Number of coupled MCMC chains
reps <- 5e2

# Step size (proposal std dev) for RWM chains
h <- 0.5

# Maximum number of coupled MCMC iterations, after lag
iter <- 1e3

# Draw starting values
sample.pi0 <- function(){rnorm(1,20,1)}

##############
# Experiment #
##############

# Lag parameter
L <- 1
# Iterations to calculate TVD bound at
bound_iters <- seq(0, 250)

runFixedLag <- function(L, bound_iters) {
  m <- 5e2
  
  # Run "reps" coupled MCMC chains at a given lag "L" 
  coupled_out <- runCoupledRWM(reps,logpi,sample.pi0,h,L,iter)
  
  ks <- seq(0, 250)
  mse_stdmcmc <- rep(NA,length(ks))
  mse_umcmc <- rep(NA,length(ks))
  
  for(i in seq_along(ks)) {
    umcmc <- allUnbiasedEstimators(coupled_out$xs,coupled_out$ys, coupled_out$taus,
                                   L, ks[i], m)
    stdmcmc <- allStandardEstimators(coupled_out$xs,
                                     ks[i], m)
    mse_umcmc[i]   <- var(umcmc)
    mse_stdmcmc[i] <- mean(stdmcmc^2)
  }
  
  plt_lim <- c(mse_umcmc, mse_stdmcmc)
  plt_lim <- c(min(plt_lim), max(plt_lim))
  # Suggested plotting code
  plot(ks, mse_umcmc, col = 'red',type = "l", log = "y",
       xlab = '"burn-in" k', 
       ylab = 'Mean-squared error', ylim = plt_lim,
       main = paste0("Lag = ", L))
  lines(ks, mse_stdmcmc, col = 'black')
  tau_quantiles <- c(quantile(coupled_out$taus, c(0.9, 0.95, 0.99)), max(coupled_out$taus))
  abline(v = tau_quantiles, col = "blue", lty = rev(seq_along(tau_quantiles)))
  
  # Plot TVD bound
  tvd_bound <- tvdCouplingBound(coupled_out$taus, L, bound_iters)
  return(tvd_bound)
}

# Lag parameter
L <- 1
# Iterations to calculate TVD bound at
bound_iters <- seq(0, 250)
tvd_bound <- runFixedLag(L,bound_iters)
plot(tvd_bound, type = "l")

########################
# Task: run the experiment for increasing values of the lag
#    1. See how the MSE behaves as the lag increases.
#    2. Compare TVD bounds as the lag increases.
########################
  

