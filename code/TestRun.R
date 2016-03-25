remove(list=ls())
library(survival)
library(gmm)
library(mvtnorm)
library(foreach)
library(doMC)

source("InfoCombine_03042016.R")
cround<-function(x){sprintf("%.2f",  round(x, 2))}


# true parameters
alpha0 = rep(0.5,4)
bx0 = c(-0.5, 0.5)

T0 = c(0.5,0.75,1,1.25)
auxdata = Sim.data(n=1, tcut0=T0) 
phi0 = auxdata$phi0


# true model 
formula0 <- as.formula("Surv(y, d) ~ z1 + z2")

# combine information
J = 1 # number of endpoints used
par0 = c(alpha0[1:J], bx0)

nrep <- 128#00 
cenp <- rep(NA, nrep) # censored percentage
fit.cox <- V.est <- matrix(NA, ncol= length(bx0), nrow=nrep) # true cox model coefficient
fit1.gmm <- fit2.gmm <-matrix(NA, ncol= length(par0), nrow=nrep)
V1.cbd<- V1.est <- matrix(NA, ncol= length(par0), nrow=nrep) # var(beta) of combined EEs and cox model

n <- 100
registerDoMC(detectCores())
tmp = foreach(k = 1:nrep, .combine=rbind) %dopar% { 
      set.seed(10101+k)
      print(k)

      tmp <- Sim.data(n=n, aux=FALSE) 
      wdata = tmp$data
      grpID = tmp$grpID
      #cenp[k] <- tmp$pcen

      
      # fit the true model
      fit <- coxph( formula0, data=wdata )
      #fit.cox[k,] <- fit$coef   
      #V.est[k,] <- diag(fit$var)

      ## combine from population survival probability
      #G = function(theta, x) getU.asym(parm=theta, x=x, T0=T0[1], phi0=phi0[,1], grpID=grpID)
      G = function(theta, x) getU.multi_asym(parm=theta, x=x, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)
       
      fit1 <- summary( gmm(g=G, x=wdata, t0 = c(alpha0[1:J],fit$coef)  ) ) # GMM estimate based on asymptotic estimating equations with cox solution as starting value.
      #fit1.gmm[k,] <- fit1$coefficient[,1]
      #V1.cbd[k,]<- fit1$coefficient[,2]^2 # the GtVG estimate of Var(beta) evaluated at the solution.

      #foo <- var(G(theta= fit1.gmm[k,], x=wdata)) # emperical covariance matrix estimate of G at the solution.
      # foo1 <- foo[1:length(fit$coef), 1:length(fit$coef)] # different blocks of th weight matrix
      # foo2 <- foo[1:length(fit$coef), 1:length(bx0)+length(fit$coef)]
      # foo3 <- foo[1:length(bx0)+length(fit$coef), 1:length(bx0)+length(fit$coef)]
      # V1.est[k,] <- diag(solve(foo1%*%solve(foo1-foo2%*%solve(foo3)%*%t(foo2))%*%t(foo1))/n) # solve for Var(beta) while avoid using derivative estimates.
      coef = c(fit$coef, fit1$coefficient[-(1:J),1])
      coef
}  

fit.cox = tmp[,1:2]
fit1.gmm = tmp[,3:4]

bias0 <- apply(fit.cox, 2, mean) - bx0
var0  <- apply(fit.cox, 2, var )
bias1 <- apply(fit1.gmm, 2, mean) - bx0
var1  <- apply(fit1.gmm, 2, var )  


sink(file = paste("out20160217-fixed-n",n,"-", option,".txt", sep="" ) )



cat("\n N =", n, "; rep =", nrep,
    "\n proportion of censoring =",  signif(mean(cenp),2),
    "\n Gamma =", cround(gamma0),
    "\n Gamma, subgroup =", cround(gamma0.SG),
    "\n True regression parameters:", bx0, bz0, 
    "\n Conventional Cox:",
    "\n\t   Bias =", cround(bias0) ,
    "\n\t     SD =", cround( sqrt(var0)) ,    
    "\n\t    ASD =", cround( sqrt(apply(V.est, 2, mean))), 
    "\n Combined information:",
    "\n\t   Bias =", cround(bias1) ,
    "\n\t     SD =", cround(sqrt(var1)) , 
    "\n\t     RE =", cround(var0/var1 ), 
    "\n\t    ASD =", cround( sqrt(apply(V1.cbd, 2, mean))),
    "\n\t Est.SD =", cround( sqrt(apply(V1.est, 2, mean))),  
    "\n Combined information (Subgroup):",
    "\n\t   Bias =", cround(bias2) ,
    "\n\t     SD =", cround(sqrt(var2)) , 
    "\n\t     RE =", cround(var0/var2), 
    "\n\t    ASD =", cround( sqrt(apply(V2.cbd, 2, mean))),
    "\n\t Est.SD =", cround( sqrt(apply(V2.est, 2, mean))) 
)

sink()

save(list=ls(), file=paste("out20160217-fixed-n",n,"-", option,".Rdata", sep="" ) )

