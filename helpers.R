# Run coupled RWM chains
runCoupledRWM <- function(reps,logpi,sample.pi0,h,L,iter) {
  xs <- matrix(NA, nrow = reps, ncol = iter + L + 1)
  ys <- matrix(NA, nrow = reps, ncol = iter + 1)
  taus <- rep(NA, reps)
  
  for(r in 1:reps) {
    x <- sample.pi0()
    y <- sample.pi0()
    out <- RWM_indepmax(x,y,h,L,iter,logpi)
    taus[r] <- out$tau
    xs[r,] <- out$xs
    ys[r,] <- out$ys
  }
  return(list("xs" = xs, "ys" = ys, "taus" = taus))
}

# Calculate the standard estimator
standardEstimator <- function(x,k,m) {
  mean(x[k:m + 1])
}

# Calculate the bias correction term
biasCorrectionTerm <- function(x,y,L,tau,k,m) {
  if (tau - 1 < k ) {
    return(0)
  } else {
    t <- k:(tau - 1)
    coeffs <- 1 + floor((t-k)/L) - ceiling(pmax(0,t - m)/L)
    div <- 1/(m - k +1)
    
    return(div * sum(coeffs * (x[t + L + 1] - y[t + 1])))
  }
}

# Calculate the unbiased estimator
unbiasedEstimator <- function(x,y,L,tau,k,m) {
  return(standardEstimator(x,k,m) + biasCorrectionTerm(x,y,L,tau,k,m))
}

# Calculate several unbiased estimators
allUnbiasedEstimators <- function(xs,ys,taus,L,k,m) {
  out <- rep(NA, nrow(xs))
  for(r in 1:nrow(xs)) {
    out[r] <- unbiasedEstimator(xs[r,],ys[r,],L,taus[r],k,m)
  }
  return(out)
}

# Calculate several standard estimators
allStandardEstimators <- function(xs,k,m) {
  out <- rep(NA, nrow(xs))
  for(r in 1:nrow(xs)) {
    out[r] <- standardEstimator(xs[r,],k,m)
  }
  return(out)
}

# Get TVD coupling bound
tvdCouplingBound <- function (tau, L, which_iters) 
{
  if (L == Inf) {
    point_estimate <- vapply(which_iters, function(t) {mean(tau > t)}, FUN.VALUE = 1)
  }
  else {
    point_estimate <- rep(NA, length(which_iters))
    i <- 0
    for (t in which_iters) {
      i <- i + 1
      point_estimate[i] <- mean(pmax(0, ceiling((tau - t)/L)))
    }
  }
  return(point_estimate)
}
