remove(list=ls())
library(methods) 
library(survival)
library(gmm)
library(mvtnorm)
library(foreach)
library(doMC)

source("InfoCombine_03042016.R")
#cround<-function(x){sprintf("%.2f",  round(x, 4))}

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1){
      J = as.numeric(args[1])
      n = as.numeric(args[2])
      nrep = as.numeric(args[3])
      maxCore = as.numeric(args[4])
} else {
      J = 1
      n = 100
      nrep = 50
      maxCore = 20
}

J1 = J+1
# true parameters
alpha0 = seq(0.4,0.7,0.1)
bx0 = c(-0.5, 0.5)

#auxdata = Sim.data(n=1, tcut0=T0, rho0=1.5) 
load("Population_data_rho1.5.Rdata")
#load("Population_data.Rdata")
T0 = c(0.5,0.75,1,1.25)
T0 = T0[c(1,3,4)]
phi0 = auxdata$phi0[,c(1,3,4)]


# true model 
formula0 <- as.formula("Surv(y, d) ~ z1 + z2")

# combine information
# J = 2 # number of time points used
par0 = c(alpha0[1:J], bx0)

#nrep <- 64*10
cenp <- rep(NA, nrep) # censored percentage
fit.cox <- V.est <- matrix(NA, ncol= length(bx0), nrow=nrep) # true cox model coefficient
fit1.gmm <- fit2.gmm <-matrix(NA, ncol= length(par0), nrow=nrep)
V1.cbd<- V1.est <- matrix(NA, ncol= length(par0), nrow=nrep) # var(beta) of combined EEs and cox model

#n <- 100
registerDoMC(min(detectCores(), maxCore))
model_fits = foreach(k = 1:nrep, .combine=rbind) %dopar% { 
      set.seed(21205+10*k)
      print(k)

      parm = c()
      for (censoring in 1:3){
            tmp <- Sim.data(n=n, aux=FALSE, sc=censoring) 
            wdata = tmp$data
            grpID = tmp$grpID
            #cenp[k] <- tmp$pcen
            
            # fit the true model
            fit <- coxph( formula0, data=wdata )

            ## combine from population survival probability
            G = function(theta, x) getU.rho.multi_asym(parm=theta, x=x, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)
            
            #fit1 <- summary( gmm(g=G, x=wdata, t0 = c(alpha0[1:J],fit$coef)  ) ) # GMM estimate based on asymptotic estimating equations with cox solution as starting value.
            fit1 = tryCatch( summary( gmm(g=G, x=wdata, t0 = c(1,alpha0[1:J],fit$coef) ) ), 
                  error = function(c) "error solving GMM." )

            if (is.character(fit1)) {
                  coef = c(fit$coef, sqrt(diag(fit$var)), rep(NA, 3*length(bx0)) )
                  message(fit1)
            } else {
                  #Sigmahat = var(G(theta= fit1$coefficient[,1], x=wdata)) # emperical covariance matrix estimate of G at the solution.
                  #Dhat = G.grad(theta = fit1$coefficient[,1], x=wdata)
                  
                  #nVarhat = tryCatch(  solve(t(Dhat)%*%solve(Sigmahat)%*%Dhat), 
                  #      error = function(c) "error inverting covariance matrix." )

                  nVarhat = "not calculating Varhat"
                  if (is.character(nVarhat)) {
                        coef = c(fit$coef, sqrt(diag(fit$var)), fit1$coefficient[-(1:J1),1], 
                              fit1$coefficient[-(1:J1),2], fit1$coefficient[1,1] )
                        #message(nVarhat)
                  } else {
                        std.err = sqrt(diag(nVarhat/n))
                        coef = c(fit$coef, sqrt(diag(fit$var)), fit1$coefficient[-(1:J),1], 
                                    fit1$coefficient[-(1:J),2], std.err[-(1:J)] )
                  }
            }
            parm = c(parm, coef)
      }
      parm
}  

fit.cox = model_fits[,1:2]
fit1.gmm = model_fits[,5:6]

bias0 <- round(apply(fit.cox, 2, mean) - bx0, 4)
var0  <- round(apply(fit.cox, 2, var ) , 4)
bias1 <- round(apply(fit1.gmm, 2, mean, na.rm=TRUE) - bx0 , 4)
var1  <- round(apply(fit1.gmm, 2, var, na.rm=TRUE) , 4) 

print(paste0("n=",n, ", J=",J, ", rho=1.5, rho_hat=TRUE" ))
print(paste("Cox: bias =",bias0[1], bias0[2], ", var = ", var0[1], var0[2], sep=" "))
print(paste("GMM: bias =",bias1[1], bias1[2], ", var = ", var1[1], var1[2], sep=" "))


save(model_fits, bx0, T0, phi0, J, file=paste("Rho1.5-Rhohat-n",n,"-J",J,".Rdata", sep="" ) )



# ------------------------------------ Bayesian GMM ------------------------------------- #
# tmp = Sim.data(n=n, aux=FALSE) 
# wdata = tmp$data
# grpID = tmp$grpID
# data.list = list(X=wdata, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)
# jump.parameter = data.frame(jump.scale=1, jump.sd1=0.02, jump.sd2=0.05)

# source("InfoCombine_03042016.R")
# mcmc.result = Bayesian.GMM(DataList=data.list, theta.init=c(alpha0[1:J],fit$coef), J=J, 
#                         nburn=1, npost=501, jump.pars=jump.parameter, 
#                         alpha.scale=5, beta.sd=10, shrinkage=TRUE)








# sink(file = paste("out20160217-fixed-n",n,"-", option,".txt", sep="" ) )

# cat("\n N =", n, "; rep =", nrep,
#     "\n proportion of censoring =",  signif(mean(cenp),2),
#     "\n Gamma =", cround(gamma0),
#     "\n Gamma, subgroup =", cround(gamma0.SG),
#     "\n True regression parameters:", bx0, bz0, 
#     "\n Conventional Cox:",
#     "\n\t   Bias =", cround(bias0) ,
#     "\n\t     SD =", cround( sqrt(var0)) ,    
#     "\n\t    ASD =", cround( sqrt(apply(V.est, 2, mean))), 
#     "\n Combined information:",
#     "\n\t   Bias =", cround(bias1) ,
#     "\n\t     SD =", cround(sqrt(var1)) , 
#     "\n\t     RE =", cround(var0/var1 ), 
#     "\n\t    ASD =", cround( sqrt(apply(V1.cbd, 2, mean))),
#     "\n\t Est.SD =", cround( sqrt(apply(V1.est, 2, mean))),  
#     "\n Combined information (Subgroup):",
#     "\n\t   Bias =", cround(bias2) ,
#     "\n\t     SD =", cround(sqrt(var2)) , 
#     "\n\t     RE =", cround(var0/var2), 
#     "\n\t    ASD =", cround( sqrt(apply(V2.cbd, 2, mean))),
#     "\n\t Est.SD =", cround( sqrt(apply(V2.est, 2, mean))) 
# )

# sink()
#save(list=ls(), file=paste("out20160217-fixed-n",n,"-", option,".Rdata", sep="" ) )

