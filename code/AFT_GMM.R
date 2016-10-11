# AFT model with subgroup survival probability info combination
library(data.table)
library(doMC)
library(foreach)
library(ggplot2)
library(gmm)
library(survival)

AFTloglik <- function(par, wdata) {
  X = as.matrix(cbind(rep(1, nrow(wdata)), wdata[, -(1:2)]))
  p = ncol(X)
  beta = par[1:p]
  sigma = par[-(1:p)]
  exb = exp(X %*% beta)
  sum(wdata$d * ((1 / sigma - 1) * log(wdata$y) -
                   log(sigma) - log(exb) / sigma) -
        (wdata$y / exb) ^ (1 / sigma))
}

AFTscore <- function(par, wdata) {
  X = as.matrix(cbind(rep(1, nrow(wdata)), wdata[, -(1:2)]))
  p = ncol(X)
  beta = par[1:p]
  sigma = par[-(1:p)]
  exb = exp(X %*% beta)
  tmp = (wdata$y / exb) ^ (1 / sigma) - wdata$d
  score.beta =  as.vector(tmp) * X / sigma
  tmp1 = wdata$d * (log(exb) - sigma - log(wdata$y))
  tmp2 = (tmp + wdata$d) * (log(wdata$y) - log(exb))
  score.sigma = (tmp1 + tmp2) / sigma ^ 2
  cbind(score.beta, score.sigma)
}

AFTSurvProb <- function(beta, sigma, x, t0) {
  exp(-(t0 * exp(-crossprod(beta, x))) ^ (1 / sigma))
}

SubGrpSurvival <- function(par, wdata, grpID, phi0, T0) {
  # Example:
  #   SubGrpSurvival(par, wdata, grpID, phi0, 0.5)
  #  
  X = as.matrix(cbind(rep(1, nrow(wdata)), wdata[, -(1:2)]))
  beta = par[-length(par)]
  sigma = par[length(par)]
  dat = cbind(grpID, X)
  surv.prb = cbind(phi0)
  
  num.grp = max(grpID)
  num.time = length(T0)
  fn = matrix(NA, nrow = nrow(X), ncol = num.grp * num.time)
  fn.index = 1
  for (t in 1:num.time) {
    for (g in 1:num.grp) {
      fn[, fn.index] = apply(dat, 1, function(x) {
        (x[1] == g) * (AFTSurvProb(beta, sigma, x[-1], T0[t]) - surv.prb[g, t])
      })
      fn.index = fn.index + 1
    }
  }
  fn
}

AFTGMMequation <- function(par, wdata, grpID, phi0, T0) {
  # Example:
  #  AFTGMMequation(par, wdata, grpID, phi0, T0 = 0.5)
  #
  score = AFTscore(par, wdata)
  aux = SubGrpSurvival(par, wdata, grpID, phi0, T0)
  cbind(score, aux)
}

FitAFTGMM <- function(wdata, grpID, phi0, T0, par.init) {
  # Example:
  #   FitAFTGMM(wdata, grpID, phi0, 0.5, par)
  #
  G = function(theta, x) {
    AFTGMMequation(par = theta, wdata = x, grpID, phi0, T0)
  }
  fit1 = tryCatch(
    summary(gmm(g = G, x = wdata, t0 = par.init, method = "BFGS")), 
    error = function(c) "error solving GMM.")
  fit1
}

SimDataAFT <- function(n.sample = 100, tcut0 = 0.5,
                       aux = FALSE, option = 1,
                       pr.cens = c(0, 0.3, 0.5, 0.7)) {
  # Example:
  #   dat.obj = SimDataAFT(100, 0.5, FALSE, 1, pr.cens = 0.3)  
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
  } else if (abs(pr.cens - 0.75) < 1e-6) {
    cens = runif(n.sample, 0, 0.836)
  } else {
    cens = Inf
    message("Choose pr.cens among: 0, 0.3, 5, 0.75.")
  }
  y = pmin(t0, cens)
  d = as.integer(t0 <= y)
  wdata = cbind(y, d, z1, z2)
  return(list(data = as.data.frame(wdata), grpID = grpid,
              pcen = 1 - sum(d) / n.sample, phi0 = phi0,
              true.par = c(beta, sigma)))
}

SimStudyAFTGMM <- function(n.sample, tcut0, n.rep, max.core,
                           option = 1, pr.cens = 0.3) {
  # Example:
  #   SimStudyAFTGMM(100, 0.5, 5, 4, 1, 0.3)
  #  
  registerDoMC(min(detectCores(), max.core))
  phi0 = SimDataAFT(n.sample, tcut0, TRUE, option, pr.cens)$phi0
  formula0 = as.formula("Surv(y, d) ~ z1 + z2")
  if (option == 2) {
    # TODO: other possible formula
  }
  # start parallel simulation
  gmm.fits = foreach(k = 1:n.rep, .combine = rbind) %dopar% { 
    set.seed(21205 + 10 * k)
    print(k)
    # get data
    dat.obj = SimDataAFT(n.sample, tcut0, FALSE, option, pr.cens)
    grpID = dat.obj$grpID
    wdata = dat.obj$data
    # fit AFT regression
    aft.obj = survreg(formula0, wdata, dist = "weibull")
    aft.mle = c(aft.obj$coefficients, aft.obj$scale)
    # fit AFT GMM
    aft.gmm = FitAFTGMM(wdata, grpID, phi0, tcut0, aft.mle)
    if (length(aft.gmm) < 5) {
      aft.gmm.coef = rep(99999, length(aft.mle))
    } else if (aft.gmm$algoInfo$convergence > 0) {
      aft.gmm.coef = rep(9999, length(aft.mle))
    } else {
      aft.gmm.coef = aft.gmm$coefficients[, 1]
    }
    output = NULL
    if (option == 1) {
      output = data.frame(method = rep(c("MLE", "GMM"), each = 4),
                          par.name = rep(c("b0", "b1", "b2", "sigma"), 2),
                          value = c(aft.mle, aft.gmm.coef))
    } else if (option == 2) {
      # TODO:
    }
    output
  }
  gmm.fits
}

# true.par = data.frame(value = c(beta, sigma),
#                       par.name = c("b0", "b1", "b2", "sigma"))

PlotSimStudy <- function(sim.study.obj, true.par) {
  dt = as.data.table(sim.study.obj)
  dt = dt[value < 99]
  g = ggplot()
  print(
    g + geom_violin(data = dt, aes(x = method, y = value),
                    draw_quantiles = c(0.025, 0.5, 0.975)) +
      facet_grid(. ~ par.name) +
      geom_hline(data = true.par,
                 aes(yintercept = value, color = "red"))
  )
}

plot(density(as.matrix(dt[par.name == "b0" & method == "MLE",
                          "value", with = FALSE])[,1]))
