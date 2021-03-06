#remove(list=ls()) 
library(survival)
library(lattice)
library(gmm)
library(matrixStats)
library(corpcor)

if (system("whoami", intern = TRUE) == "dengdetian0603") {
  source('~/Documents/JHSPH/Research/CY.Huang/Code/SurvivalGMM/code/SurvivalDataSimulation.R')
} else {
  source('/users/ddeng/ThesisTopic/SurvivalGMM/code/SurvivalDataSimulation.R')
}


# single survival endpoint, multiple subgroups
getU.asym <- function(parm, x=wdata, T0, phi0, grpID)
{
      gdata <- x; N <- nrow(gdata)
      # score function
      alpha = parm[1]
      beta = parm[-1]
      gx <- as.matrix( gdata[,-1*(1:2)] )
      ebx <- exp( gx %*% beta )
      S0 <- sapply( gdata$y, function(u, t=gdata$y, a=ebx){ sum(a*(t>=u)) }) # n*s^(0)( t=y_i, beta) for i=1,...,N
      #baseline hazard
      dLam <- gdata$d/S0 #assume no ties; this needs to be taken care of later
      ### dN(y_i)/nS^(0)(y_i)  ####### use dN(y_i)/n for E[dN(t)]|y_i !!!!!

      findInt <- function(u, tt=gdata$y, dh ){ sum(dh[tt<=u]) } ## \int_{t <= u}^ dh(t)

      S1 <- matrix(NA, ncol=ncol(gx), nrow=N)
      foo1 <- NULL
      for(kk in 1: ncol(gx)){
            S1[,kk] <- sapply( gdata$y, function(u, t=gdata$y, a=gx[,kk]*ebx){ sum(a*(t>=u)) }) ## the kk-th element of n*S^(1)( Y_i, beta) for i=1,...,N
            tmp <- gdata$d *(gx[,kk] - S1[,kk]/S0 )   ## \int_{t} [Z_i(t=y_i) - s^(1)/s^(0)]dN_i(t)
            tmp <- tmp - (gx[,kk] * ebx * sapply(gdata$y, findInt, dh=dLam) - 
                  ebx* sapply(gdata$y, findInt, dh=dLam * S1[,kk]/S0) ) 
            ## 
            foo1 <- cbind(foo1, tmp)
      } 
      U1<- foo1


      # accumulative hazard
      tmp <- U2<- NULL
      tmp = sapply(gdata$y, min, T0)
      U2 = N*dLam * (gdata$y<=T0) - N*ebx * sapply(tmp, findInt, dh=dLam/S0 ) + sapply(rep(T0, N), findInt, dh=dLam) - alpha

      # auxiliary info
      foo1 <-  tmp <- U3 <- NULL
      grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
      for (j in 1:length(phi0))
      {
            tmp =  (grpID == grpval[j])*(exp(-alpha*ebx) - phi0[j])
            foo1 = cbind(foo1, tmp)
      }
      U3 = foo1

      return(cbind(U1, U2, U3))
}


# gradient of U bar
getGrad = function(parm, x, T0, phi0, grpID)
{
      index = order(x$y, decreasing=TRUE)
      gdata = x[index, ]
      grpID = as.vector(grpID)[index]
      
      N = nrow(gdata)
      J = length(T0)
      if (J > 1) {
            K = nrow(phi0)
            dU2_dalpha = diag(rep(-1, J))
      } else {
            K = length(phi0) 
            phi0 = matrix(phi0, ncol=1)
            dU2_dalpha = -1
      }

      alpha = parm[1:J]
      beta = parm[-(1:J)]
      p = length(beta)

      gx = as.matrix( gdata[,-1*(1:2)] )
      ebx = exp( gx %*% beta )
      #S0 <- sapply( gdata$y, function(u, t=gdata$y, a=ebx){ sum(a*(t>=u)) }) # n*s^(0)( t=y_i, beta) for i=1,...,N
      xebx = diag(ebx[,1]) %*% gx
      xxt =  t(apply(gx,1,function(x) {x%*%t(x)})) # each xxt matrix is compressed as a row vector
      xxtebx = diag(ebx[,1]) %*% xxt

      nS0 = cumsum(ebx)
      nS1 = colCumsums(xebx)
      nS2 = colCumsums(xxtebx)

      #baseline hazard
      dLam = gdata$d/nS0 #assume no ties; this needs to be taken care of later
      ### dN(y_i)/nS^(0)(y_i)  ####### use dN(y_i)/n for E[dN(t)]|y_i !!!!!

      findInt <- function(u, tt=gdata$y, dh ){ sum(dh[tt<=u]) } ## \int_{t <= u}^ dh(t)  

      dU1_dalpha = matrix(0, nrow=p, ncol=J)
      

      tmp = NULL
      dU3_dalpha = matrix(0, nrow=J*K, ncol=J)
      dU3_dbeta = matrix(0, nrow=J*K, ncol=p)

      grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
      for (j in 1:J){     
            for (k in 1:K){
                  tmp =  -(grpID == grpval[k])*(exp(-alpha[j]*ebx)*ebx)
                  dU3_dalpha[(j-1)*K+k,j] = mean(tmp) 
                  tmp = -(grpID == grpval[k])*(exp(-alpha[j]*ebx)*alpha[j]*ebx)
                  #print(tmp)
                  tmp = diag(tmp[,1]) %*% gx
                  dU3_dbeta[(j-1)*K+k,] = apply(tmp, 2, mean)
            }
      }

      ##  dU1_dbeta = 
      S1S1t =  t(apply(nS1/N,1,function(x) {x%*%t(x)}))
      tmp1 = diag(dLam)%*%(nS2/N - diag(N/nS0)%*%S1S1t)
      du1.part1 = matrix(apply(tmp1, 2, sum), nrow=p, ncol=p)

      du1.part2 = matrix(0, nrow=p, ncol=p)
      for (i in 1:N){
            for (ii in i:N){
                  tmp = -xxt[i,] + (nS2[ii,]+2*as.vector(nS1[ii,]%*%t(gx[i,])))/nS0[ii] - 2*S1S1t[ii,]*N*N/nS0[ii]^2
                  du1.part2 = du1.part2 + matrix(ebx[i]*dLam[ii]*tmp, nrow=p, ncol=p)
            }
      }
      du1.part2 = du1.part2/N
      dU1_dbeta = -du1.part1 - du1.part2

      ##  TO be corrected: dU2_dbeta = 
      du2.part1 = t(sapply(T0, function(t){ apply(diag((gdata$y<=t)/nS0*dLam)%*%nS1, 2, sum) } )) * (-2)

      lessthanTj = sapply(T0, function(t) as.numeric(gdata$y<=t))
      #
      du2.part2 = matrix(NA, nrow=J, ncol=p)
      du2.part3 = matrix(NA, nrow=J, ncol=p)
      for (j in 1:J){
            tmp = dLam*lessthanTj[,j]/nS0
            tmp1 = diag(tmp*N/nS0) %*% nS1
            tmp2 = diag(ebx[,1]) %*% colCumsums(tmp1[N:1,])[N:1,]
            tmp3 = diag(cumsum(tmp[N:1])[N:1]) %*% xebx
            du2.part2[j,] = apply(tmp3, 2, sum)
            du2.part3[j,] = apply(tmp2, 2, mean)*2
      }
      dU2_dbeta = du2.part1 - du2.part2 + du2.part3

      Gradient = matrix(NA, nrow=p+J+J*K, ncol=J+p)
      Gradient[1:p, 1:J] = dU1_dalpha
      Gradient[1:p, (J+1):(J+p)] = dU1_dbeta
      Gradient[(p+1):(p+J), 1:J] = dU2_dalpha
      Gradient[(p+1):(p+J), (J+1):(J+p)] = dU2_dbeta
      Gradient[(p+J+1):(p+J+J*K), 1:J] = dU3_dalpha
      Gradient[(p+J+1):(p+J+J*K), (J+1):(J+p)] = dU3_dbeta
      return(Gradient)
}

# Gbar = function(theta, x) apply(getU.multi_asym(parm=theta, x=x, T0=T0, phi0=phi0, grpID=grpID), 2, mean)
# getGrad(parm=c(0.3,0.6,-0.5,0.1), wdata, T0, phi0, grpID) -> aa
# theta = c(0.3,0.6,-0.5,0.1)
# attr( numericDeriv(quote(Gbar(theta,wdata)), c("theta")) , "gradient") -> aaa

# multiple survival endpoint, multiple subgroups
getU.multi_asym <- function(parm, x=wdata, T0, phi0, grpID) # phi be a J by K matrix
{
      gdata <- x; N <- nrow(gdata)
      J = length(T0)
      if (J > 1) K = nrow(phi0)
      else {
            K = length(phi0)
            phi0 = matrix(phi0, ncol=1)
      }
      # score function
      alpha = parm[1:J]
      beta = parm[-(1:J)]
      gx <- as.matrix( gdata[,-1*(1:2)] )
      ebx <- exp( gx %*% beta )
      S0 <- sapply( gdata$y, function(u, t=gdata$y, a=ebx){ sum(a*(t>=u)) }) # n*s^(0)( t=y_i, beta) for i=1,...,N
      #baseline hazard
      dLam <- gdata$d/S0 #assume no ties; this needs to be taken care of later
      ### dN(y_i)/nS^(0)(y_i)  ####### use dN(y_i)/n for E[dN(t)]|y_i !!!!!

      findInt <- function(u, tt=gdata$y, dh ){ sum(dh[tt<=u]) } ## \int_{t <= u}^ dh(t)

      S1 <- matrix(NA, ncol=ncol(gx), nrow=N)
      U1 = matrix(NA, ncol=ncol(gx), nrow=N)
      for(kk in 1: ncol(gx)){
            S1[,kk] <- sapply( gdata$y, function(u, t=gdata$y, a=gx[,kk]*ebx){ sum(a*(t>=u)) }) ## the kk-th element of n*S^(1)( Y_i, beta) for i=1,...,N
            tmp <- gdata$d *(gx[,kk] - S1[,kk]/S0 )   ## \int_{t} [Z_i(t=y_i) - s^(1)/s^(0)]dN_i(t)
            U1[, kk] <- tmp - (gx[,kk] * ebx * sapply(gdata$y, findInt, dh=dLam) - 
                  ebx* sapply(gdata$y, findInt, dh=dLam * S1[,kk]/S0) ) 
      } 

      # accumulative hazard # check the correctness ########################
      tmp <- NULL
      U2 = matrix(NA, ncol=J, nrow=N)
      for(j in 1:J){
            tmp = sapply(gdata$y, min, T0[j])
            U2[,j] = N*dLam * (gdata$y<=T0[j]) - N*ebx * sapply(tmp, findInt, dh=dLam/S0 ) + sapply(rep(T0[j], N), findInt, dh=dLam) - alpha[j]
      }
      
      # auxiliary info
      tmp <- NULL
      # grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
      grpval = 1:K
      U3 = matrix(NA, ncol = J*K, nrow=N)
      c = 1
      for (j in 1:J){
            for (k in 1:K){
                  tmp =  (grpID == grpval[k])*(exp(-alpha[j]*ebx) - phi0[k,j])
                  U3[, c] = tmp
                  c = c + 1
            }
      }

      return(cbind(U1, U2, U3))
}
# getU.multi_asym(parm=c(alpha0[1:J],fit$coef), x=wdata, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)



# multiple survival endpoint, multiple subgroups, with allowing rho != 1
getU.rho.multi_asym <- function(parm, x=wdata, T0, phi0, grpID) # phi be a J by K matrix
{
      gdata <- x; N <- nrow(gdata)
      J = length(T0)
      if (J > 1) K = nrow(phi0)
      else {
            K = length(phi0)
            phi0 = matrix(phi0, ncol=1)
      }
      # score function
      J1 = J+1
      rho = parm[1]
      alpha = parm[2:J1]
      beta = parm[-(1:J1)]

      gx <- as.matrix( gdata[,-1*(1:2)] )
      ebx <- exp( gx %*% beta )
      S0 <- sapply( gdata$y, function(u, t=gdata$y, a=ebx){ sum(a*(t>=u)) }) # n*s^(0)( t=y_i, beta) for i=1,...,N
      #baseline hazard
      dLam <- gdata$d/S0 #assume no ties; this needs to be taken care of later
      ### dN(y_i)/nS^(0)(y_i)  ####### use dN(y_i)/n for E[dN(t)]|y_i !!!!!

      findInt <- function(u, tt=gdata$y, dh ){ sum(dh[tt<=u]) } ## \int_{t <= u}^ dh(t)

      S1 <- matrix(NA, ncol=ncol(gx), nrow=N)
      U1 = matrix(NA, ncol=ncol(gx), nrow=N)
      for(kk in 1: ncol(gx)){
            S1[,kk] <- sapply( gdata$y, function(u, t=gdata$y, a=gx[,kk]*ebx){ sum(a*(t>=u)) }) ## the kk-th element of n*S^(1)( Y_i, beta) for i=1,...,N
            tmp <- gdata$d *(gx[,kk] - S1[,kk]/S0 )   ## \int_{t} [Z_i(t=y_i) - s^(1)/s^(0)]dN_i(t)
            U1[, kk] <- tmp - (gx[,kk] * ebx * sapply(gdata$y, findInt, dh=dLam) - 
                  ebx* sapply(gdata$y, findInt, dh=dLam * S1[,kk]/S0) ) 
      } 

      # accumulative hazard # check the correctness ########################
      tmp <- NULL
      U2 = matrix(NA, ncol=J, nrow=N)
      for(j in 1:J){
            tmp = sapply(gdata$y, min, T0[j])
            U2[,j] = N*dLam * (gdata$y<=T0[j]) - N*ebx * sapply(tmp, findInt, dh=dLam/S0 ) + sapply(rep(T0[j], N), findInt, dh=dLam) - alpha[j]
      }
      
      # auxiliary info
      tmp <- NULL
      # grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
      grpval = 1:K
      U3 = matrix(NA, ncol = J*K, nrow=N)
      c = 1
      for (j in 1:J){
            for (k in 1:K){
                  tmp =  (grpID == grpval[k])*(exp(-rho*alpha[j]*ebx) - phi0[k,j])
                  U3[, c] = tmp
                  c = c + 1
            }
      }

      return(cbind(U1, U2, U3))
}
# getU.multi_asym(parm=c(alpha0[1:J],fit$coef), x=wdata, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)

