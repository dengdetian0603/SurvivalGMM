
SimSurvival1 <- function(n.sample = 100, tcut0 = 0.5,
                         aux = FALSE, option = 1,
                         pr.cens = c(0, 0.3, 0.5, 0.75)) {
  # Simulate survival time from weibull distribution, satisfy both Cox and AFT
  # without interaction term 
  # Example:
  #   dat.obj = SimDataAFT(500, 0.5, FALSE, 1, pr.cens = 0.3)  
  #
  Xcut <- function(w, v){ 
    cbind((w <= qnorm(0.5) & v == 0), (w > qnorm(0.5) & v == 0))  
  }
  beta = c(0.2, 0.5, -0.5)
  sigma = 0.5
  phi0 = NULL
  if (option == 2) {
    # TODO:
  }
  if (aux) {
    n <- 500000
    z1 <- rnorm( n )
    z2 <- rbinom( n, 1, 0.5 )
    ZB = beta[1] + z1 * beta[2] + z2 * beta[3]
    Acut0 <- Xcut(z1, z2 )
    phi0 = matrix(NA, ncol = length(tcut0), nrow = ncol(Acut0))
    t0 = sapply(exp(ZB), rweibull, n = 1, shape = 1/sigma)
    for (k in 1:ncol(Acut0)) {
      for (j in 1:length(tcut0)) {
        phi0[k, j] = sapply(tcut0[j], function(u, t = t0[Acut0[, k]]) {
          mean(t > u) 
        })
      }
    }
  }
  z1 = rnorm(n.sample)
  z2 = rbinom(n.sample, 1, 0.5)
  Acut0 <- Xcut(z1, z2 )
  grpid = Acut0  %*% c(1, 2)
  ZB = beta[1] + z1 * beta[2] + z2 * beta[3]
  t0 = sapply(exp(ZB), rweibull, n = 1, shape = 1/sigma)
  if (abs(pr.cens - 0) < 1e-6) {
    cens = Inf
  } else if (abs(pr.cens - 0.3) < 1e-6) {
    cens = runif(n.sample, 0, 3.2)
  } else if (abs(pr.cens - 0.5) < 1e-6) {
    cens = runif(n.sample, 0, 1.725)
  } else if (abs(pr.cens - 0.45) < 1e-6) {
    cens = abs(rnorm(n.sample, .4, .5)) * (rbinom(n.sample, 1, 0.5) * 4 + 1)
  } else if (abs(pr.cens - 0.75) < 1e-6) {
    cens = runif(n.sample, 0, 0.836)
  } else if (abs(pr.cens - 0.8) < 1e-6) {
    cens = rexp(n.sample, 3)
  } else {
    cens = Inf
    message("Choose pr.cens among: 0, 0.3, 0.5, 0.75.")
  }
  y = pmin(t0, cens)
  d = as.integer(t0 <= y)
  wdata = cbind(y, d, z1, z2)
  return(list(data = as.data.frame(wdata), grpID = grpid,
              pcen = 1 - sum(d) / n.sample, phi0 = phi0,
              true.par = c(beta, sigma)))
}


SimSurvival2 <- function(n.sample = 100, tcut0 = 0.5,
                         aux = FALSE, option = 1, aux.rho = 1,
                         pr.cens = c(0, 0.3, 0.5, 0.75)) {
  # Simulate survival time from weibull distribution, satisfy both Cox and AFT
  # with interaction term 
  # Example:
  #   dat.obj = SimDataAFT(500, 0.5, FALSE, 1, pr.cens = 0.3)  
  #
  Xcut <- function(w, v){ 
    cbind((w <= qnorm(0.5) & v == 0), (w > qnorm(0.5) & v == 0))  
  }
  beta = c(0.2, 0.5, -0.5, 0.5)
  sigma = 0.5
  phi0 = NULL
  if (option == 2) {
    # TODO:
  }
  if (aux) {
    n <- 500000
    z1 <- rnorm( n )
    z2 <- rbinom( n, 1, 0.5 )
    ZB = beta[1] + z1 * beta[2] + z2 * beta[3] + z1 * z2 * beta[4]
    Acut0 <- Xcut(z1, z2 )
    phi0 = matrix(NA, ncol = length(tcut0), nrow = ncol(Acut0))
    t0 = sapply(exp(ZB)/sqrt(aux.rho), rweibull, n = 1, shape = 1/sigma)
    for (k in 1:ncol(Acut0)) {
      for (j in 1:length(tcut0)) {
        phi0[k, j] = sapply(tcut0[j], function(u, t = t0[Acut0[, k]]) {
          mean(t > u) 
        })
      }
    }
  }
  z1 = rnorm(n.sample)
  z2 = rbinom(n.sample, 1, 0.5)
  Acut0 <- Xcut(z1, z2 )
  grpid = Acut0  %*% c(1, 2)
  ZB = beta[1] + z1 * beta[2] + z2 * beta[3] + z1 * z2 * beta[4]
  t0 = sapply(exp(ZB), rweibull, n = 1, shape = 1/sigma)
  if (abs(pr.cens - 0) < 1e-6) {
    cens = Inf
  } else if (abs(pr.cens - 0.3) < 1e-6) {
    cens = runif(n.sample, 0, 3.2)
  } else if (abs(pr.cens - 0.5) < 1e-6) {
    cens = runif(n.sample, 0, 1.725)
  } else if (abs(pr.cens - 0.45) < 1e-6) {
    cens = abs(rnorm(n.sample, .4, .5)) * (rbinom(n.sample, 1, 0.5) * 4 + 1)
  } else if (abs(pr.cens - 0.75) < 1e-6) {
    cens = runif(n.sample, 0, 0.836)
  } else if (abs(pr.cens - 0.8) < 1e-6) {
    cens = rexp(n.sample, 3)
  } else {
    cens = Inf
    message("Choose pr.cens among: 0, 0.3, 0.5, 0.75.")
  }
  y = pmin(t0, cens)
  d = as.integer(t0 <= y)
  wdata = cbind(y, d, z1, z2)
  return(list(data = as.data.frame(wdata), grpID = grpid,
              pcen = 1 - sum(d) / n.sample, phi0 = phi0,
              true.par = c(beta, sigma)))
}

SimSurvival3 <- function(n.sample = 100, tcut0 = 0.5,
                         aux = FALSE, option = 1, aux.rho = 1,
                         pr.cens = c(0, 0.3, 0.5, 0.75)) {
  # Simulate survival time from weibull distribution, satisfy both Cox and AFT
  # with interaction term 
  # Example:
  #   dat.obj = SimDataAFT(500, 0.5, FALSE, 1, pr.cens = 0.3)  
  #
  Xcut <- function(w, v){ 
    cbind((w <= qnorm(0.5) & v == 0), (w > qnorm(0.5) & v == 0))  
  }
  beta = c(0.2, 0.5, -0.5, 0) # no interaction
  sigma = 0.5
  phi0 = NULL
  if (option == 2) {
    # TODO:
  }
  if (aux) {
    n <- 500000
    z1 <- rnorm( n )
    z2 <- rbinom( n, 1, 0.5 )
    ZB = beta[1] + z1 * beta[2] + z2 * beta[3] + z1 * z2 * beta[4]
    Acut0 <- Xcut(z1, z2 )
    phi0 = matrix(NA, ncol = length(tcut0), nrow = ncol(Acut0))
    t0 = sapply(exp(ZB)/sqrt(aux.rho), rweibull, n = 1, shape = 1/sigma)
    for (k in 1:ncol(Acut0)) {
      for (j in 1:length(tcut0)) {
        phi0[k, j] = sapply(tcut0[j], function(u, t = t0[Acut0[, k]]) {
          mean(t > u) 
        })
      }
    }
  }
  z1 = rnorm(n.sample)
  z2 = rbinom(n.sample, 1, 0.5)
  Acut0 <- Xcut(z1, z2 )
  grpid = Acut0  %*% c(1, 2)
  ZB = beta[1] + z1 * beta[2] + z2 * beta[3] + z1 * z2 * beta[4]
  t0 = sapply(exp(ZB), rweibull, n = 1, shape = 1/sigma)
  if (abs(pr.cens - 0) < 1e-6) {
    cens = Inf
  } else if (abs(pr.cens - 0.3) < 1e-6) {
    cens = runif(n.sample, 0, 3.2)        # mean(cens <= t0)
  } else if (abs(pr.cens - 0.5) < 1e-6) {
    cens = runif(n.sample, 0, 1.725)
  } else if (abs(pr.cens - 0.45) < 1e-6) {
    cens = abs(rnorm(n.sample, .4, .5)) * (rbinom(n.sample, 1, 0.5) * 4 + 1)
  } else if (abs(pr.cens - 0.75) < 1e-6) {
    cens = runif(n.sample, 0, 0.836)
  } else if (abs(pr.cens - 0.8) < 1e-6) {
    cens = rexp(n.sample, 3)
  } else {
    cens = Inf
    message("Choose pr.cens among: 0, 0.3, 0.5, 0.75.")
  }
  y = pmin(t0, cens)
  d = as.integer(t0 <= y)
  wdata = cbind(y, d, z1, z2)
  return(list(data = as.data.frame(wdata), grpID = grpid,
              pcen = 1 - sum(d) / n.sample, phi0 = phi0,
              true.par = c(beta, sigma)))
}




