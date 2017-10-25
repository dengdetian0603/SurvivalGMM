library(aftgee)
library(dplyr)

if (system("whoami", intern = TRUE) == "dengdetian0603") {
  source('~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/AFT_GMM_logrank.R')
} else {
  source("/users/ddeng/ThesisTopic/SurvivalGMM/code/AFT_GMM_logrank.R")
}
# Rcpp::sourceCpp("~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/AFT_GMM_logrank.cpp")

t.star = c(0.5, 1, 1.2)
dat.obj0 = SimDataAFT(150, t.star, TRUE, 1, pr.cens = 0.45)
surv.prob = dat.obj0$phi0


dat.obj = SimDataAFT(800, t.star, FALSE, 1, pr.cens = 0.45)  

# beta.vec = c(0.5, -0.5)
dat.mat = as.matrix(dat.obj$data)
grp.id = dat.obj$grpID[, 1]
grp.id0 = grp.id
# ------------------------------------------------------------------------------
surv.moments = survMoments4(dat.mat, grp.id0, surv.prob, t.star)
colMeans(surv.moments)

Moments = FullMoments(beta.vec, dat.mat, grp.id0, surv.prob, t.star)
# Moments = FullMoments2(beta.vec, dat.mat, grp.id, surv.prob, t.star)
# cor(Moments)
colMeans(Moments)

round(cor(Moments, surv.moments), 4)

plot(Moments[, 3], Moments[, 7], col = dat.obj$grpID + 1)

gmmFunc <- function(theta) {
  ymat = FullMoments(beta.vec = theta, dat.mat = dat.mat,
                     grp.id, surv.prob, t.star)
  colMeans(ymat)
}

z1 = seq(0.45, 0.55, by = 5 * 1e-4)
y = sapply(z1, function(x) gmmFunc(c(x, -0.5)))
plot(z1, y[1, ], type = "l")
lines(z1, y[2, ], col = 2)
lines(z1, y[3, ], col = 3)
lines(z1, y[4, ], col = 4)
lines(z1, y[5, ], col = 5)

D = jacobian(gmmFunc, x = beta.vec, method = "Richardson",
             method.args = list(eps = 1e-4, d = 0.0001,
                              zero.tol = .Machine$double.eps * 1e+5,
                              r = 6, v = 2, show.details = TRUE))
precMat = t(D) %*% solve(Omega) %*% D
solve(precMat)
# ------------------------------------------------------------------------------
# use FullMoments for subgroup equation constructed by non-censored data
# use FullMoments2 for subgroup equation constructed by censored data
gmmEq <- function(theta, x) {
	FullMoments(beta.vec = theta, dat.mat = x, grp.id, surv.prob, t.star)
}

gmmEq2 <- function(theta, x) {
  a = FullMoments(beta.vec = theta, dat.mat = x, grp.id, surv.prob, t.star)
  cbind(a, surv.moments)
}



fit0 = aftsrr(Surv(y, d) ~ z1 + z2, data = dat.obj$data,
              rankWeights = "gehan", method = "sm",
              B = 50, variance = c("ISMB", "ZLMB"))
summary(fit0)

gfit1 = suppressWarnings(gmm(g = GehanMoments, x = as.matrix(dat.obj$data),
                             t0 = coef(fit0), method = "BFGS", type = "two"))
gfit1

t0 = proc.time()
gfit2 = suppressWarnings(gmm(g = gmmEq2, x = as.matrix(dat.obj$data),
                             t0 = coef(fit0), method = "BFGS", type = "two"))
proc.time() - t0
gfit2

rbind(coef(fit0), coef(gfit1), coef(gfit2))

# ----------------------------------------------------------------------------------
t.star = c(0.5, 0.9, 1.2)
dat.obj0 = SimDataAFT(150, t.star, TRUE, 1, pr.cens = 0.45)
surv.prob = dat.obj0$phi0

res = SimStudyAFTGMM(n.sample = 120, tcut0 = t.star, n.rep = 150,
                     max.core = 4, option = 1, pr.cens = 0.45,
                     surv.prob = surv.prob, cover.prob = FALSE,
                     IPCW.moments = TRUE)

res %>% filter(abs(value) < 2) %>% group_by(method, par.name) %>%
  do(data.frame(mean = mean(.$value), sd = sd(.$value)))

true.par = data.frame(value = c(0.5, -0.5), par.name = c("b1", "b2"))
PlotSimStudy(res, true.par)

# ----------------------------------------------------------------------------------
W = as.matrix(gfit2$w0)
var.beta = VarBeta(beta.vec = coef(gfit2), dat.mat, grp.id, surv.prob, t.star,
                   B = 1000, W)
std.err.beta = sqrt(diag(var.beta))
std.err.beta

colMeans(res[res$method == "GMM" & res$par.name == "b1",
         c("std.err", "cover.by.ci")])
colMeans(res[res$method == "GMM" & res$par.name == "b2",
         c("std.err", "cover.by.ci")])
