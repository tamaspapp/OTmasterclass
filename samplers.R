#############################################
# RWM
#############################################
RWM <- function (x,h,iter,logpi) {
  d <- length(x)
  if(d!=1){
    stop("Fix dimension to 1 for proper output handling.")
  }
  logpi_x <- logpi(x)

  accepts <- 0
  tau <- -1
  
  xs <- rep(NA, iter+1)
  xs[1] <- x
  
  for(i in 1:iter) {
    # Generate X-proposal
    xdot <- rnorm(d)
    xp <- x + h * xdot
    logpi_xp <- logpi(xp)
    
    u_acc.log <- log(runif(1))
    if (u_acc.log < logpi_xp - logpi_x) {
      x <- xp
      logpi_x <- logpi_xp
    }
    xs[i+1] <- x
  }
  return(list("x" = x, "xs" = xs))
}



################################################
# RWM: Maximal coupling w/ independent residuals
################################################
RWM_indepmax <- function (x,y,
                          h,
                          L,    # Lag
                          iter, # Max number of iterations AFTER lag
                          logpi) {
  d <- length(x)
  if(d!=1){
    stop("Fix dimension to 1 for proper output handling.")
  }
  
  xs <- rep(NA, 1 + L + iter)
  ys <- rep(NA, 1 + iter)
  
  init_out <- RWM(x,h,L,logpi)
  
  x <- init_out$x
  xs[0:L + 1] <- init_out$xs
  ys[1] <- y
  
  logpi_x <- logpi(x)
  logpi_y <- logpi(y)
  
  accepts <- 0
  tau <- -1
  
  for(i in 1:iter) {
    # Generate X-proposal
    z <- (x-y)/h
    xdot <- rnorm(d)
    xp <- x + h * xdot
    logpi_xp <- logpi(xp)
    
    prob_cplprop.log <- 0.5*(sum(xdot^2) - sum((xdot+z)^2)) # chance of coupling the proposals
    u_cplprop.log <- log(runif(1))
    if(u_cplprop.log < prob_cplprop.log){ # Propose Yp = Xp
      yp <- xp
      logpi_yp <- logpi_xp
    } else { # Rejection-sample Yp from non-overlap
      reject <- T
      while(reject) {
        ydot <- rnorm(d)
        if(log(runif(1)) > -0.5*sum(z^2)  + sum(ydot * z) ) {
          reject <- F
        }
      }
      yp <- y + h * ydot
      logpi_yp <- logpi(yp)
    }
    u_acc.log <- log(runif(1))
    if (u_acc.log < logpi_xp - logpi_x) {
      x <- xp
      logpi_x <- logpi_xp
    }
    if (u_acc.log < logpi_yp - logpi_y) {
      y <- yp
      logpi_y <- logpi_yp
    }
    xs[i + L+1] <- x
    ys[i + 1]   <- y
    
    if(identical(x,y)) {
      tau <- i
      break
    }
  }
  
  if(tau == -1) {stop("Chains not coupled. Increase iteration count.")}
  if(tau < iter) { # Do final (iter - tau) iterations
    final_out <- RWM(x,h,iter-tau,logpi)
    xs[(tau+1):iter + L + 1] <- final_out$xs[-1]
  }
  
  return(list("xs" = xs,
              "ys" = ys,
              "tau" = tau))
}



#############################################
# Maximal coupling w/ independent residuals
#############################################