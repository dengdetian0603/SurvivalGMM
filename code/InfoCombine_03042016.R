
remove(list=ls()) 
library(survival)
library(lattice)
library(gmm)

Sim.data = function(b1=-0.5, b2=0.5, option=2, sc=2, rho0=1, n=100, aux = TRUE, tcut0=0.5){
      phi0 = NULL
      Acut0 = NULL

      if (option ==1 & sc ==2 ){ ucen <- 2.6721; cenp  <- 30 }#censoring pr=30%
      if (option ==1 & sc ==3 ){ ucen <- 1.2253; cenp  <- 50 }#censoring pr=50%
      if (option ==2 & sc ==2 ){ ucen <- 2.6995; cenp  <-30} #censoring pr=30%
      if (option ==2 & sc ==3 ){ ucen <- 1.5647; cenp  <-50} #censoring pr=30% 

      #=====================================#  
      #   generate auxiliary information    #
      #=====================================#  
      Xcut <- function(w, v){ 
            cbind((w<=qnorm(0.5) & v==0), (w >qnorm(0.5) & v==0) )  
      }

      if (aux){
            n <- 500000
            z1 <- rnorm( n )
            z2 <- rbinom( n, 1, 0.5 )  

            if (option==1){
                  # lambda_0(t)=1 #Lam0(t)=(t/scale)^shape
                  sfoo <- exp(-1*(log(rho0)+b1*z1 + b2*z2))
                  t0 <- sapply(sfoo, rweibull, n=1, shape=1)  
            }
            if (option==2){ 
                  # lambda_0(t)=2t
                  sfoo <- exp(-1*(log(rho0)+b1*z1 + b2*z2)/2)
                  t0 <- sapply(sfoo, rweibull, n=1, shape=2)      
            }

            Acut0 <- Xcut(z1, z2 )
            # S( tcut[j] | Acut[,k] )   
            phi0 = matrix(NA, nrow=length(tcut0), ncol=ncol(Acut0))
            for( k in 1: ncol(Acut0) ){
                  for (j in 1:length(tcut0)){
                        phi0[k,j] = sapply(tcut0[j], function(u, t=t0[Acut0[,k]]){ mean(t>u) } )
                  }
            } 
      }

      #=====================================#  
      #         simulate data set           #
      #=====================================#  
      z1 <- rnorm( n )
      z2 <- rbinom( n, 1, 0.5 )  

      #Lam0(t)=(t/scale)^shape
      if (option==1){
            # lambda_0(t)=1
            sfoo <- exp(-1*(b1*z1 + b2*z2))
            t0 <- sapply(sfoo, rweibull, n=1, shape=1)  
      }
      if (option==2){ 
            # lambda_0(t)=2t
            sfoo <- exp(-1*(b1*z1 + b2*z2)/2)
            t0 <- sapply(sfoo, rweibull, n=1, shape=2)      
      }
      if (option==3){ 
            bz <- b1*z1 + b2*z2
            t0<-NULL
            for(mm in 1:n){
                  FF <- function(x, kk) { 1-exp(-(x-2)^3/3/2*exp(kk)) }
                  uu <- runif(1,0,1)
                  FFF <- function (x, u=uu, k) {FF(x, kk=k)-u}    
                  T0 <-  uniroot(FFF, c(0,100), k=bz[mm])$root
                  t0 <- c(t0, T0)
            }
      }

      if (sc==1) {cen=10^10} else {
            cen <- runif(n, 0, ucen)
      }

      y <- apply( cbind(t0, cen), 1, min)
      d <- as.integer( t0<=y )
      wdata <- data.frame( cbind(y=y, d=d, z1=z1, z2=z2) ) 
      pcen <- 1- mean(d)

      tcut0 <- 0.5
      Acut0 <- Xcut(z1, z2 )
      grpid = Acut0%*%c(1,2)

      return(list(data = wdata, grpID= grpid, pcen = pcen, phi0 = phi0))
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
      U2 = dLam * (gdata$y<=T0) - ebx * sapply(tmp, findInt, dh=dLam/S0 ) + sapply(rep(T0, N), findInt, dh=dLam) - alpha

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
getGrad = function(parm, x, U, T0, phi0, grpID)
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

      dU1_dalpha = matrix(0, nrow=length(beta), ncol=1)
      dU2_dalpha = matrix(1, nrow=1, ncol=1)

      K = length(phi0)
      tmp = NULL
      dU3_dalpha = matrix(NA, nrow=K, ncol=1)
      dU3_dbeta = matrix(NA, nrow=K, ncol=length(beta))
      grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
      for (k in 1:K)
      {
            tmp =  -(grpID == grpval[k])*(exp(-alpha*ebx)*ebx)
            dU3_dalpha[k,1] = mean(tmp) 
            tmp = -(grpID == grpval[k])*(exp(-alpha*ebx)*alpha*ebx)
            tmp = diag(tmp) %*% gx
            dU3_dbeta[k,] = apply(tmp, 2, mean)
      }

## TODO:
      dU1_dbeta = 

}



# multiple survival endpoint, multiple subgroups
getU.multi_asym <- function(parm, x=wdata, T0, phi0, grpID) # phi be a J by K matrix
{
      gdata <- x; N <- nrow(gdata)
      J = length(T0)
      K = nrow(phi0)

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

      # accumulative hazard
      tmp <- NULL
      U2 = matrix(NA, ncol=J, nrow=N)
      for(j in 1:J){
            tmp = sapply(gdata$y, min, T0[j])
            U2[,j] = dLam * (gdata$y<=T0[j]) - ebx * sapply(tmp, findInt, dh=dLam/S0 ) + sapply(rep(T0[j], N), findInt, dh=dLam) - alpha[j]
      }
      
      # auxiliary info
      tmp <- NULL
      grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
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





#### Combining auxilliary information project
#*# confirm GMM method works 
## extend to multiple survival end-points
## try Bayesian GMM with 
## combine median follow-up time

#### Imputation Censoring project
## Efron's re-dristribution to the right (EM) principle: self-consistency algorithm.
## 


