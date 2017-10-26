# simulation study
# CoxPH, AFT induced GMM estimator
# with or without interaction
# 1 or 2 landmark times  
library(aftgee)
library(data.table)
library(doMC)
library(foreach)
library(gmm)
library(lattice)
library(matrixStats)
library(methods)
library(survival)



dir = ifelse(system("whoami", intern = TRUE) == "dengdetian0603",
             '~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/',
             '/users/ddeng/ThesisTopic/SurvivalGMM/') 

source(paste0(dir, 'code/SurvivalDataSimulation.R'))
source(paste0(dir, 'code/CoxPH_GMM.R'))
source(paste0(dir, 'code/AFT_GMM_logrank.R'))

#------------------------------------------------------------------------------
setwd(dir)
if (file.exists('./code/population_data.Rdata')) {
  load('./code/population_data.Rdata')
} else {
  t.star = c(0.5, 1, 1.2)
  dat.obj0.a = SimSurvival1(150, t.star, TRUE, 1, pr.cens = 0.3)
  dat.obj0.b = SimSurvival2(150, t.star, TRUE, 1, pr.cens = 0.3)
  save(t.star, dat.obj0.a, dat.obj0.b, file = './code/population_data.Rdata')
}

tid = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(tid)) tid = 0

# =============================================================================
grid = expand.grid(num.T = 1:3, pr.cens = c(0.3, 0.5, 0.75))

## no interaction
surv.prob = dat.obj0.a$phi0
formula0 = as.formula("Surv(y, d) ~ z1 + z2")

out = c()
for (i in 1:nrow(grid)) {
  num.T = grid$num.T[i]
  t.star.i = t.star[1:num.T]
  surv.prob.i = cbind(surv.prob[, 1:num.T])

  dat.obj = SimSurvival1(150, t.star.i, FALSE, 1,
                         pr.cens = grid$pr.cens[i])
  dat.mat = as.matrix(dat.obj$data)
  grp.id = dat.obj$grpID[, 1]
  grp.id0 = grp.id    
  
  # AFT Gehan ~~~~~~~~~~~~~~~~~~~~~~~~~~
  gmmEq <- function(theta, x) {
    FullMoments(beta.vec = theta, dat.mat = x, grp.id, surv.prob.i, t.star.i)
  }

  fit0 = aftsrr(formula0, data = dat.obj$data,
              rankWeights = "gehan", method = "sm",
              B = 50, variance = c("ZLMB"))

  gfit2 = suppressWarnings(gmm(g = gmmEq, x = dat.mat,
                               t0 = coef(fit0), method = "BFGS", type = "two"))

  W = as.matrix(gfit2$w0)
  var.beta = VarBeta(beta.vec = coef(gfit2), dat.mat, grp.id,
                     surv.prob.i, t.star.i, B = 1000, W)
  std.err.beta = sqrt(diag(var.beta))
        
  res.1 = c(tid,'AFT', num.T, grid$pr.cens[i],
            coef(fit0), coef(gfit2), std.err.beta, NA, NA)

  # Cox PH ~~~~~~~~~~~~~~~~~~~~~~~~~~
  fit = coxph(formula0, data = dat.obj$data)

  G = function(theta, x) {
    getU.multi_asym(parm = theta, x = x, T0 = t.star.i,
                    phi0 = surv.prob.i, grpID = grp.id)
  }
  G.grad = function(theta, x) {
    getGrad(parm = theta, x = x, T0 = t.star.i,
            phi0 = surv.prob.i, grpID = grp.id) 
  }

  alpha0 = seq(0.4,0.7,0.1) # initial values of nuisance parameter
  fit1 = tryCatch(gmm(g = G, x = as.data.frame(dat.mat),
                      t0 = c(alpha0[1:num.T], fit$coef)),
                  error = function(c) "error solving Cox-GMM.")

  if (is.character(fit1)) {
    val = c(fit$coef, rep(NA, 6))
    message(fit1)
  } else {
    Sigmahat = var(G(theta = fit1$coefficient, x = as.data.frame(dat.mat)))
    Dhat = G.grad(theta = fit1$coefficient, x = as.data.frame(dat.mat))
    
    nVarhat = tryCatch(solve(t(Dhat) %*% solve(Sigmahat) %*% Dhat), 
                       error = function(c) "error inverting cov matrix." )
    if (is.character(nVarhat)) {
          val = c(fit$coef,
                  fit1$coefficient[-(1:num.T)], 
                  NA, NA,
                  sqrt(diag(fit1$vcov))[-(1:num.T)])
          message(nVarhat)
    } else {
          std.err = sqrt(diag(nVarhat/nrow(dat.mat)))[-(1:num.T)]
          val = c(fit$coef,
                  fit1$coefficient[-(1:num.T)], 
                  std.err,
                  sqrt(diag(fit1$vcov))[-(1:num.T)])          
    }
  }
  res.2 = c(tid, 'Cox', num.T, grid$pr.cens[i], val)

  out = rbind(out, res.1, res.2)
}

colnames(out) = c('tid', 'Model', 'Num_T', 'Pr_Censor',
                  'NoAux_Beta1_est', 'NoAux_Beta2_est',
                  'Aux_Beta1_est', 'Aux_Beta2_est',
                  'Aux_Beta1_stderr', 'Aux_Beta2_stderr',
                  'Aux_Beta1_stderr_', 'Aux_Beta2_stderr_')

# --------------------------------------------------------------------------
file.name1 = "./experiments/NoInteraction.csv"

if (file.exists(file.name1)) {
  write.table(out, sep = ",", row.names = FALSE,
              col.names = FALSE,
              file = file.name1, append = TRUE)
} else {
  write.table(out, sep = ",", row.names = FALSE,
              col.names = TRUE,
              file = file.name1)
}





