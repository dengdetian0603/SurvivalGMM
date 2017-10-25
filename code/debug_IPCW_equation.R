library(aftgee)
library(dplyr)

source('~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/AFT_GMM_logrank.R')
Rcpp::sourceCpp("~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/AFT_GMM_logrank.cpp")

t.star = c(0.5, 1)
dat.obj0 = SimDataAFT(150, t.star, TRUE, 1, pr.cens = 0.45)
surv.prob = dat.obj0$phi0


dat.obj = SimDataAFT(800, t.star, FALSE, 1, pr.cens = 0.45)  

# beta.vec = c(0.5, -0.5)
dat.mat = as.matrix(dat.obj$data)
grp.id = dat.obj$grpID[, 1]
grp.id0 = grp.id

# ------------------------------------------------------------------------------

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

colMeans(Psi.mat)
# ------------------------------------------------------------------------------
k = 1; t = 1
mean( d.vec * (grp.id == k) * exp(gamma.mat[, k])) / mean(grp.id == k)

mean(d.vec * (y.vec >= t.star[t] & grp.id == k) * exp(gamma.mat[, k]) -
     surv.prob[k, t] * (grp.id == k))

mean(surv.prob[k, t] * (grp.id == k))

mean(d.vec * (y.vec >= t.star[t] & grp.id == k) * exp(gamma.mat[, k])) / mean(grp.id == k)

plot(y.vec, exp(-gamma.mat[, 1]), type = "l")
lines(y.vec, exp(-gamma.mat[, 2]), col = 2)

cbind(d.vec, y.vec, grp.id, exp(-gamma.mat[, 1]))[y.vec >= t.star[1] & grp.id == 1, ] -> b1
cbind(d.vec, y.vec, grp.id, exp(-gamma.mat[, 2]))[y.vec >= t.star[1] & grp.id == 2, ] -> b2
plot(b1[, 2], b1[, 1]/b1[, 4], col = b1[, 1] + 1)
plot(b2[, 2], b2[, 1]/b2[, 4], col = b2[, 1] + 1)

surv.prob.hat = matrix(NA, nrow = num.grp, ncol = num.time)
for (k in 1:num.grp) {
  for (t in 1:num.time) {
    surv.prob.hat[k, t] = 
      mean(d.vec * (y.vec >= t.star[t] & grp.id == k) *
           exp(gamma.mat[, k])) / mean(grp.id == k)
    
  }
}
surv.prob.hat
surv.prob

# ------------------------------------------------------------------------------
beta.vec = c(0.5, -0.5)
Moments = FullMoments(beta.vec, dat.mat, grp.id0, surv.prob, t.star)

round(cor(cbind(Moments, Psi.mat)), 4)
