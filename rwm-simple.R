source("samplers.R")
source("helpers.R")

#########
# Setup #
#########

# Target log-density
logpi <- function(x){-0.5*sum(x^2)}

# Number of coupled MCMC chains
reps <- 1e3

# Step size (proposal std dev) for RWM chains
h <- 0.5

# Maximum number of coupled MCMC iterations, after lag
iter <- 1e3

# Draw starting values
sample.pi0 <- function(){rnorm(1,20,1)}

##############
# Experiment #
##############

# Run "reps" coupled MCMC chains at lag "L = 1" 
coupled_out <- runCoupledRWM(reps,logpi,sample.pi0,h,L = 1,iter)

# Calculate standard and unbiased MCMC estimators for given "k", "m"
# NOTES: 
# - "m" needs to be less than "L + iter"
# - "k" can be at most "m"

m <- 5e2
k <- 1.7e2

umcmc <- allUnbiasedEstimators(coupled_out$xs,coupled_out$ys, coupled_out$taus,
                               L=1, k, m)
stdmcmc <- allStandardEstimators(coupled_out$xs,
                                 k, m)
mse_stdmcmc <- mean(stdmcmc^2)
mse_umcmc   <- var(umcmc)

mse_stdmcmc
mse_umcmc

########################################################################################
# Task: fix "m" to 500 or 1000, say, and calculate MSE for k = 0 to largest meeting time
########################################################################################

m <- 5e2

## Suggested plotting code
# plot(ks, mse_umcmc, col = 'red',type = "l", log = "y",
#      xlab = '"burn-in" k', 
#      ylab = 'MSE (logarithmic scale)')
# lines(ks, mse_stdmcmc, col = 'black')
