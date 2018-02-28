# AFT model with subgroup survival probability info combination
library(data.table)
library(doMC)
library(foreach)
library(ggplot2)
library(gmm)
library(survival)

if (system("whoami", intern = TRUE) == "dengdetian0603") {
  Rcpp::sourceCpp("~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/AFT_GMM_logrank.cpp")
  source('~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/SurvivalDataSimulation.R')
} else {
  Rcpp::sourceCpp("/users/ddeng/ThesisTopic/SurvivalGMM/code/AFT_GMM_logrank.cpp")
  source('/users/ddeng/ThesisTopic/SurvivalGMM/code/SurvivalDataSimulation.R')
}

GehanMoments <- function(beta.vec, dat.mat) {
# gval = GehanMoments(c(0.5, -0.5), as.matrix(dat.obj$data))
# colMeans(gval)
# gfit = gmm(g = GehanMoments, x = as.matrix(dat.obj$data),
#           t0 = c(0.5, -0.5), method = "BFGS", type = "two")
#
  X.mat = cbind(dat.mat[, -(1:2)])
  N = nrow(dat.mat)

  e.vec = log(dat.mat[, 1]) - X.mat %*% cbind(beta.vec)
  e.order = order(e.vec)
  order0 = (1:N)[e.order]
  e.vec = e.vec[e.order]
  X.mat = X.mat[e.order, ]
  d.vec = dat.mat[e.order, 2]
  
  Xd.mat = X.mat * d.vec
  XdR.mat = Xd.mat * (N : 1)

  G1.mat = t(apply(Xd.mat, 2, cumsum)/N)
  G2.vec = colSums(XdR.mat)/(N^2)
  G3.mat = t(XdR.mat/N)
  G4.mat = t(d.vec/N * apply(X.mat[N:1, ], 2, cumsum)[N:1, ])
  G5.mat = t(cumsum(d.vec)/N * X.mat)
  G6.vec = rowMeans(G4.mat)

  G.mat = t(G1.mat - G2.vec + G3.mat - G4.mat - G5.mat + G6.vec)
  G.mat[order0, ]
}


survMoments0 <- function(rho = NULL, beta.vec, dat.mat,
                         grp.id, surv.prob, t.star) {
  # survival info condition through breslow estimator
  # allowing flexible baseline hazard 
  X.mat = cbind(dat.mat[, -(1:2)])
  N = nrow(dat.mat)
  Xbeta = X.mat %*% cbind(beta.vec)
  
  e.vec = log(dat.mat[, 1]) - Xbeta
  e.order = order(e.vec)
  order0 = (1:N)[e.order]
  e.vec = e.vec[e.order]
  X.mat = X.mat[e.order, ]
  d.vec = dat.mat[e.order, 2]
  
  # Sub-group survivial probabilities
  # ----------------------------------------------------------------------------
  surv.prob = cbind(surv.prob)
  num.grp = max(grp.id)
  num.time = length(t.star)
  grp.id = as.vector(grp.id)[e.order]
  
  epsilon.mat = sapply(log(t.star),  function(x) x - Xbeta[e.order, ])
  dR.vec = d.vec / (N:1)
  alpha.mat = hazardMat(cumsum(dR.vec), e.vec, epsilon.mat)
  ndR2.vec = N * d.vec / ((N:1)^2)
  
  if (length(rho) < 1) {
    Psi.mat = psiMat(cumsum(ndR2.vec), e.vec, d.vec,
                     epsilon.mat, alpha.mat, grp.id, surv.prob)    
  } else {
    Psi.mat = psiMatFlex(rho, cumsum(ndR2.vec), e.vec, d.vec,
                         epsilon.mat, alpha.mat, grp.id, surv.prob)
  }
  Psi.mat[order0, ]
}


survMoments2 <- function(dat.mat, grp.id0, surv.prob, t.star) {
  N = nrow(dat.mat)
  y.vec0 = dat.mat[, 1]
  y.order = order(y.vec0)
  order0 = (1:N)[y.order]
  y.vec = y.vec0[y.order]
  dc.vec = 1 - dat.mat[y.order, 2]
  grp.id = as.vector(grp.id0)[y.order]
  
  surv.prob = cbind(surv.prob)
  num.grp = max(grp.id)
  num.time = length(t.star)
  
  Gamma.mat = matrix(NA, nrow = N, ncol = num.grp * num.time)
  for (k in 1:num.grp) {
    pr.grp.k = mean(grp.id == k)
    y.vec.k = y.vec[grp.id == k]
    dc.vec.k = dc.vec[grp.id == k]
    Sc.k = c(length(y.vec.k):1, 1e-16)/N
    Sc.u.k = sapply(y.vec, function(x) {
      Sc.k[binSearch2(y.vec.k, x) + 1]
    })
    for (t in 1:num.time) {
      Sc.tk = Sc.k[binSearch2(y.vec.k, t.star[t]) + 1]
      Nc.tk = dc.vec.k * (y.vec.k <= t.star[t])
      gamma.tk = sum(Nc.tk/Sc.k[1:length(y.vec.k)])/N # cencoring cumulative hazard
      
      common.factor = Sc.tk * exp(gamma.tk)/pr.grp.k
      summand.a = (y.vec >= t.star[t] & grp.id == k)/Sc.tk -
                  (grp.id == k)/pr.grp.k
      summand.b = (y.vec <= t.star[t] & grp.id == k) * dc.vec/Sc.u.k
      
      tmp1 = outer(1:N, 1:N, FUN = ">=") * (grp.id == k)
      tmp2 = (y.vec <= t.star[t] & grp.id == k) * dc.vec/(Sc.u.k^2)
      summand.c = tmp1 %*% cbind(tmp2)/N
      
      Gamma.mat[, (k - 1)*num.time + t] =
        as.vector(common.factor * (1 + summand.a + summand.b - summand.c[, 1]) 
                  - surv.prob[k, t])
    }
  }
  Gamma.mat[order0, ]
}

survMoments3 <- function(dat.mat, grp.id0, surv.prob, t.star) {
  # subgroup survival info quations derived from IPCW estimator
  N = nrow(dat.mat)
  y.vec0 = dat.mat[, 1]
  y.order = order(y.vec0)
  order0 = (1:N)[y.order]
  y.vec = y.vec0[y.order]
  d.vec = dat.mat[y.order, 2]
  dc.vec = 1 - d.vec
  grp.id = as.vector(grp.id0)[y.order]
  
  surv.prob = cbind(surv.prob)
  num.grp = max(grp.id0)
  num.time = length(t.star)
  
  # cummulative hazard of censored data, C.
  gamma.mat = matrix(NA, nrow = N, ncol = num.grp)
  Psi.mat = matrix(NA, nrow = N, ncol = num.grp * num.time)
  for (k in 1:num.grp) {
    y.vec.k = y.vec[grp.id == k]
    dc.vec.k = dc.vec[grp.id == k]
    Sc.k = c(length(y.vec.k):1, 1e-16)/N
    Sc.u.k = sapply(y.vec, function(x) {
      Sc.k[binSearch2(y.vec.k, x) + 1]
    })
    
    V.mat.k = matrix(NA, nrow = N, ncol = N)
    for (i in 1:N) {
      Nc.ik = dc.vec.k * (y.vec.k <= y.vec[i])
      gamma.mat[i, k] = sum(Nc.ik/Sc.k[1:length(y.vec.k)])/N
      
      tmp0 = (y.vec <= y.vec[i] & grp.id == k)
      summand.b = tmp0 * dc.vec/Sc.u.k
      tmp1 = outer(1:N, 1:N, FUN = ">=") * (grp.id == k)
      tmp2 = tmp0 * dc.vec/(Sc.u.k^2)
      summand.c = tmp1 %*% cbind(tmp2)/N
      
      V.mat.k[i, ] = gamma.mat[i, k] + as.vector(summand.b - summand.c)
    }
    
    exp.gamma.k = exp(gamma.mat[, k])
    for (t in 1:num.time) {
      A.vec.tk = d.vec * (y.vec >= t.star[t] & grp.id == k)
      h0 = A.vec.tk * exp.gamma.k * V.mat.k
      h.mat = h0 + t(h0)
      hii.vec = diag(h.mat)
      EU = (sum(h.mat) - sum(hii.vec))/(N^2 - N)
        
      Psi.mat[, (k - 1)*num.time + t] =
        A.vec.tk * (1 - gamma.mat[, k]) * exp.gamma.k +
        rowMeans(h.mat) - 0.5 * EU + 0.5/N * hii.vec -
        surv.prob[k, t] * (grp.id == k)
    }
  }
  Psi.mat[order0, ]
}

survMoments4 <- function(dat.mat, grp.id0, surv.prob, t.star) {
  # subgroup survival info quations derived from IPCW estimator
  # assume independent censoring
  N = nrow(dat.mat)
  y.vec0 = dat.mat[, 1]
  y.order = order(y.vec0)
  order0 = (1:N)[y.order]
  y.vec = y.vec0[y.order]
  d.vec = dat.mat[y.order, 2]
  dc.vec = 1 - d.vec
  grp.id = as.vector(grp.id0)[y.order]
  
  surv.prob = cbind(surv.prob)
  num.grp = max(grp.id0)
  num.time = length(t.star)
  
  # cummulative hazard of censored data, C.
  gamma.vec = rep(NA, N)
  V.mat.k = matrix(NA, nrow = N, ncol = N)
  Sc.u = (N:1)/N
  for (i in 1:N) {
    Nc.i = dc.vec * (y.vec <= y.vec[i])
    gamma.vec[i] = sum(Nc.i/(N:1))
    
    tmp0 = (y.vec <= y.vec[i])
    summand.b = tmp0 * dc.vec/Sc.u
    tmp1 = outer(1:N, 1:N, FUN = ">=")
    tmp2 = tmp0 * dc.vec/(Sc.u^2)
    summand.c = tmp1 %*% cbind(tmp2)/N
    
    V.mat.k[i, ] = gamma.vec[i] + as.vector(summand.b - summand.c)
  }
  exp.gamma.k = exp(gamma.vec)
  
  Psi.mat = matrix(NA, nrow = N, ncol = num.grp * num.time)
  for (k in 1:num.grp) {
    y.vec.k = y.vec[grp.id == k]
    dc.vec.k = dc.vec[grp.id == k]
    
    for (t in 1:num.time) {
      A.vec.tk = d.vec * (y.vec >= t.star[t] & grp.id == k)
      h0 = A.vec.tk * exp.gamma.k * V.mat.k
      h.mat = h0 + t(h0)
      hii.vec = diag(h.mat)
      EU = (sum(h.mat) - sum(hii.vec))/(N^2 - N)
      
      Psi.mat[, (k - 1)*num.time + t] =
        A.vec.tk * (1 - gamma.vec) * exp.gamma.k +
        rowMeans(h.mat) - 0.5 * EU + 0.5/N * hii.vec -
        surv.prob[k, t] * (grp.id == k)
    }
  }
  Psi.mat[order0, ]
}

# ------------------------------------------------------------------------------

FullMoments <- function(beta.vec, dat.mat, grp.id, surv.prob, t.star) {
  X.mat = cbind(dat.mat[, -(1:2)])
  N = nrow(dat.mat)
  Xbeta = X.mat %*% cbind(beta.vec)

  e.vec = log(dat.mat[, 1]) - Xbeta
  e.order = order(e.vec)
  order0 = (1:N)[e.order]
  e.vec = e.vec[e.order]
  X.mat = X.mat[e.order, ]
  d.vec = dat.mat[e.order, 2]
  
  Xd.mat = X.mat * d.vec
  XdR.mat = Xd.mat * (N:1)
  # Gehan Moments
  # ----------------------------------------------------------------------------
  G1.mat = t(apply(Xd.mat, 2, cumsum)/N)
  G2.vec = colSums(XdR.mat)/(N^2)
  G3.mat = t(XdR.mat/N)
  G4.mat = t(d.vec/N * apply(X.mat[N:1, ], 2, cumsum)[N:1, ])
  G5.mat = t(cumsum(d.vec)/N * X.mat)
  G6.vec = rowMeans(G4.mat)

  G.mat = t(G1.mat - G2.vec + G3.mat - G4.mat - G5.mat + G6.vec)
  # Sub-group survivial probabilities
  # ----------------------------------------------------------------------------
  surv.prob = cbind(surv.prob)
  num.grp = max(grp.id)
  num.time = length(t.star)
  grp.id = as.vector(grp.id)[e.order]

  epsilon.mat = sapply(log(t.star),  function(x) x - Xbeta[e.order, ])
  dR.vec = d.vec / (N:1)
  alpha.mat = hazardMat(cumsum(dR.vec), e.vec, epsilon.mat)
  ndR2.vec = N * d.vec / ((N:1)^2)

  Psi.mat = psiMat(cumsum(ndR2.vec), e.vec, d.vec, epsilon.mat, alpha.mat, 
                   grp.id, surv.prob)

  cbind(G.mat, Psi.mat)[order0, ]
}

FullMomentsFlex <- function(rho, beta.vec, dat.mat, grp.id, surv.prob, t.star) {
  # allowing flexible baseline hazard 
  X.mat = cbind(dat.mat[, -(1:2)])
  N = nrow(dat.mat)
  Xbeta = X.mat %*% cbind(beta.vec)
  
  e.vec = log(dat.mat[, 1]) - Xbeta
  e.order = order(e.vec)
  order0 = (1:N)[e.order]
  e.vec = e.vec[e.order]
  X.mat = X.mat[e.order, ]
  d.vec = dat.mat[e.order, 2]
  
  Xd.mat = X.mat * d.vec
  XdR.mat = Xd.mat * (N:1)
  # Gehan Moments
  # ----------------------------------------------------------------------------
  G1.mat = t(apply(Xd.mat, 2, cumsum)/N)
  G2.vec = colSums(XdR.mat)/(N^2)
  G3.mat = t(XdR.mat/N)
  G4.mat = t(d.vec/N * apply(X.mat[N:1, ], 2, cumsum)[N:1, ])
  G5.mat = t(cumsum(d.vec)/N * X.mat)
  G6.vec = rowMeans(G4.mat)
  
  G.mat = t(G1.mat - G2.vec + G3.mat - G4.mat - G5.mat + G6.vec)
  # Sub-group survivial probabilities
  # ----------------------------------------------------------------------------
  surv.prob = cbind(surv.prob)
  num.grp = max(grp.id)
  num.time = length(t.star)
  grp.id = as.vector(grp.id)[e.order]
  
  epsilon.mat = sapply(log(t.star),  function(x) x - Xbeta[e.order, ])
  dR.vec = d.vec / (N:1)
  alpha.mat = hazardMat(cumsum(dR.vec), e.vec, epsilon.mat)
  ndR2.vec = N * d.vec / ((N:1)^2)
  
  Psi.mat = psiMatFlex(rho, cumsum(ndR2.vec), e.vec, d.vec,
                       epsilon.mat, alpha.mat, grp.id, surv.prob)
  
  cbind(G.mat, Psi.mat)[order0, ]
}

# ==============================================================================

SimDataAFT <- SimSurvival1

SimStudyAFTGMM <- function(n.sample, tcut0, n.rep, max.core,
                           option = 1, pr.cens = 0.3, surv.prob = NULL,
                           cover.prob = FALSE, IPCW.moments = FALSE) {
  # Example:
  #   SimStudyAFTGMM(100, 0.5, 5, 4, 1, 0.3)
  #  
  registerDoMC(min(detectCores(), max.core))
  if (length(surv.prob) < 1) {
    surv.prob = SimDataAFT(n.sample, tcut0, TRUE, option, pr.cens)$phi0
  }
  formula0 = as.formula("Surv(y, d) ~ z1 + z2")
  if (option == 2) {
    # TODO: other possible formula
  }
  # start parallel simulation
  gmm.fits = foreach(k = 1:n.rep, .combine = rbind) %dopar% { 
    set.seed(21205 - 7 * k)
    print(k)
    # get data
    dat.obj = SimDataAFT(n.sample, tcut0, FALSE, option, pr.cens)
    grpID = dat.obj$grpID

    # fit AFT log-rank smoothed regression
    aft.srr = aftsrr(formula0, data = dat.obj$data,
                     rankWeights = "gehan", method = "sm",
                     B = 0, variance = c("ISMB", "ZLMB"))
    aft.spmle = coef(aft.srr)
    # fit AFT GMM
    if (IPCW.moments) {
      ipcw.moments = survMoments4(dat.obj$data, grpID, surv.prob, tcut0)
      gmmEq <- function(theta, x) {
        a = GehanMoments(beta.vec = theta, dat.mat = x)
        cbind(a, ipcw.moments)
      }
    } else {
      gmmEq <- function(theta, x) {
        FullMoments(beta.vec = theta, dat.mat = x,
                    grp.id = dat.obj$grpID, surv.prob = surv.prob,
                    t.star = tcut0)
      }
    }
    aft.gmm = tryCatch(suppressWarnings(
      gmm(g = gmmEq, x = as.matrix(dat.obj$data), t0 = aft.spmle,
          method = "BFGS", type = "cue")),
      error = function(c) "error solving GMM.")
    if (length(aft.gmm) < 5) {
      aft.gmm.coef = rep(99999, length(aft.spmle))
    } else if (aft.gmm$algoInfo$convergence > 0) {
      aft.gmm.coef = rep(9999, length(aft.spmle))
    } else {
      aft.gmm.coef = coef(aft.gmm)
    }
    std.err = rep(NA, 4)
    covered = rep(NA, 4)
    if (cover.prob & length(aft.gmm) > 5) {
      W = as.matrix(aft.gmm$w0)
      var.beta = VarBeta(beta.vec = aft.gmm.coef, as.matrix(dat.obj$data),
                         dat.obj$grpID, surv.prob, tcut0, B = 900, W)
      std.err.beta = sqrt(diag(var.beta))
      std.err[3:4] = std.err.beta
      lower = aft.gmm.coef - 1.96 * std.err.beta
      upper = aft.gmm.coef + 1.96 * std.err.beta
      covered[3] = lower[1] <= 0.5 & upper[1] >= 0.5
      covered[4] = lower[2] <= -0.5 & upper[2] >= -0.5
    }
    output = NULL
    if (option == 1) {
      output = data.frame(method = rep(c("MLE", "GMM"), each = 2),
                          par.name = rep(c("b1", "b2"), 2),
                          value = c(aft.spmle, aft.gmm.coef),
                          std.err = std.err, cover.by.ci = covered)
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

# plot(density(as.matrix(dt[par.name == "b0" & method == "MLE",
#                           "value", with = FALSE])[,1]))


GradMoments <- function(beta.vec, dat.mat,
                        grp.id, surv.prob, t.star, B = 20) {
  N = nrow(dat.mat)
  p = length(beta.vec)
  DZsum = 0
  gVal = colSums(FullMoments(beta.vec, dat.mat,
                             grp.id, surv.prob, t.star))/sqrt(N)
  for (i in 1:B) {
    Z = rnorm(p, 0, 1)
    gPerturb = colSums(FullMoments(beta.vec + Z/sqrt(N), dat.mat,
                                   grp.id, surv.prob, t.star))/sqrt(N)
    DZ = cbind(gPerturb - gVal) %*% rbind(Z)
    DZsum = DZsum + DZ
  }
  DZsum/B
}


VarBeta <- function(beta.vec, dat.mat, grp.id, surv.prob,
                    t.star, B = 20, W = NULL) {
  # sieveGMM with non-smooth moments.
  g = FullMoments(beta.vec, dat.mat,
                  grp.id, surv.prob, t.star)
  V = kernHAC(lm(g ~ 1), sandwich = FALSE)
  dG = GradMoments(beta.vec, dat.mat,
                   grp.id, surv.prob, t.star, B)
  if (length(W) < 1) {
    varBeta = solve(t(dG) %*% V %*% dG)/nrow(dat.mat)
  } else {
    invGWG = solve(t(dG) %*% W %*% dG)
    GWVWG = t(dG) %*% W %*% V %*% W %*% dG
    varBeta = invGWG %*% GWVWG %*% invGWG/nrow(dat.mat)
  }
  varBeta
}


GradMomentsRho <- function(theta, dat.mat,
                           grp.id, surv.prob, t.star, B = 20) {
  N = nrow(dat.mat)
  p = length(theta)
  DZsum = 0
  gVal = colSums(FullMomentsFlex(theta[1], theta[-1], dat.mat,
                                 grp.id, surv.prob, t.star))/sqrt(N)
  for (i in 1:B) {
    Z = rnorm(p, 0, 1)
    theta.per = theta + Z/sqrt(N)
    gPerturb = colSums(FullMomentsFlex(theta.per[1], theta.per[-1], dat.mat,
                                       grp.id, surv.prob, t.star))/sqrt(N)
    DZ = cbind(gPerturb - gVal) %*% rbind(Z)
    DZsum = DZsum + DZ
  }
  DZsum/B
}


VarBetaRho <- function(theta, dat.mat, grp.id, surv.prob,
                       t.star, B = 20, W = NULL) {
  # sieveGMM with non-smooth moments.
  g = FullMomentsFlex(theta[1], theta[-1], dat.mat,
                      grp.id, surv.prob, t.star)
  V = kernHAC(lm(g ~ 1), sandwich = FALSE)
  dG = GradMomentsRho(theta, dat.mat,
                      grp.id, surv.prob, t.star, B)
  if (length(W) < 1) {
    varBeta = solve(t(dG) %*% V %*% dG)/nrow(dat.mat)
  } else {
    invGWG = solve(t(dG) %*% W %*% dG)
    GWVWG = t(dG) %*% W %*% V %*% W %*% dG
    varBeta = invGWG %*% GWVWG %*% invGWG/nrow(dat.mat)
  }
  varBeta
}







