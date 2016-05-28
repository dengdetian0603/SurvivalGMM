#remove(list=ls()) 
library(survival)
library(lattice)
library(gmm)
library(matrixStats)
library(corpcor)

Sim.data = function(b1=-0.5, b2=0.5, option=2, sc=2, rho0=1, n=100, aux = TRUE, tcut0=0.5){
      phi0 = NULL
      Acut0 = NULL

      if (option ==1 & sc ==2 ){ ucen <- 2.6721; cenp  <- 30 }#censoring pr=30%
      if (option ==1 & sc ==3 ){ ucen <- 1.2253; cenp  <- 50 }#censoring pr=50%
      if (option ==2 & sc ==2 ){ ucen <- 2.6995; cenp  <-30} #censoring pr=30%
      if (option ==2 & sc ==3 ){ ucen <- 1.5647; cenp  <-50} #censoring pr=50% 

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
            phi0 = matrix(NA, ncol=length(tcut0), nrow=ncol(Acut0))
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
      grpval = sort(unique(grpID))[-1] # grpID == 0 means it does not fit in any subgroup
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





#### Combining auxilliary information project
#*# confirm GMM method works 
#*# extend to multiple survival end-points
#*# improve optimization stablity with exact gradient
## allow constant difference between sample and population accumulative hazard
## try Bayesian GMM
## combine median follow-up time

#### Imputation Censoring project
## Efron's re-dristribution to the right (EM) principle: self-consistency algorithm.
## Compare with Jon's Censoring Unbiased Survival Tree method.

#------------------- Bayesian GMM ------------------------------#
GMM.likelihood = function(theta, X, T0, phi0, grpID, shrinkage=TRUE, logscale=TRUE){
      N = nrow(X)
      G = getU.multi_asym(parm=theta, x=X, T0=T0, phi0=phi0, grpID=grpID)
      Gbar = apply(G, 2, mean)
      if (shrinkage) {
            invSigma = invcov.shrink(G, verbose=FALSE)
            Qn = Gbar%*%invSigma%*%Gbar*N
      } else {
            Sigma = var(G)
            Qn = Gbar%*%solve(Sigma, Gbar)*N
      }
      Qn = as.vector(Qn)
      if (logscale){
            return(-0.5*Qn)
      } else {
            return(exp(-0.5*Qn))
      }
}
# GMM.likelihood(theta = fit1$coefficients[,1], X=wdata, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)

Get.prior = function(theta, J, alpha.shape, alpha.scale, beta.sd, logscale=TRUE){
      # gamma process prior for cumulative hazard function      
      alpha.increment = (c(theta[1:J],0)-c(0,theta[1:J]))[1:J]
      prior.alpha = dgamma(x=alpha.increment, shape=alpha.shape, scale=alpha.scale, log=logscale)
      prior.beta = dnorm(x=theta[-(1:J)], mean=0, sd = beta.sd, log=logscale)

      if (logscale){
            return(sum(prior.alpha) + sum(prior.beta))
      } else {
            return(prod(prior.alpha)*prod(prior.beta))
      }
}
# Get.prior(theta = fit1$coefficients[,1], J=J, alpha.shape, alpha.scale=2, beta.sd=10)

Jump2Next = function(theta0, J, iter, alphalist, jump.sd){
      alpha.new = alphalist[iter, ]
      beta.new =  theta0[-(1:J)] + rnorm(length(theta0)-J, 0, jump.sd)
      return( c(alpha.new, beta.new) )
}

Jump2Next.b = function(theta0, J, jump.sd1, jump.sd2){
      alpha.increment = (c(theta0[1:J],0)-c(0,theta0[1:J]))[1:J]
      new.increment = abs(alpha.increment + rnorm(J, 0, jump.sd1))
      alpha.new = cumsum(new.increment)

      beta.new =  theta0[-(1:J)] + rnorm(length(theta0)-J, 0, jump.sd2)
      return( c(alpha.new, beta.new) )
}

# DataList = list(X=wdata, T0=T0[1:J], phi0=phi0[,1:J], grpID=grpID)
Bayesian.GMM = function(DataList, theta.init, J, nburn, npost, jump.pars, alpha.scale, beta.sd, shrinkage=TRUE){
      nchain = npost + nburn
      acceptance = rep(0, nchain)
      
      jump.scale = jump.pars$jump.scale
      jump.sd1 = jump.pars$jump.sd1
      jump.sd2 = jump.pars$jump.sd2

      post.chain = matrix(NA, nrow=nchain, ncol=length(theta.init))
      post.chain[1,] = theta.init

      surv.prob = DataList$phi0
      if (J<2) {
            surv.prob = matrix(surv.prob)
      }
      R = apply(-log(surv.prob), 2, mean)/jump.scale
      jump.shape = (c(R,0)-c(0,R))[1:J]
      R = apply(-log(surv.prob), 2, mean)/alpha.scale
      Alpha.shape = (c(R,0)-c(0,R))[1:J]

      log.density = function(parm){
            #print(parm)
            GMM.likelihood(theta=parm, X=DataList$X, T0=DataList$T0, 
                  phi0=DataList$phi0, grpID=DataList$grpID, shrinkage=shrinkage, logscale=TRUE) + 
            Get.prior(theta=parm, J=J, alpha.shape=Alpha.shape, alpha.scale=alpha.scale, beta.sd=beta.sd, logscale=TRUE)
      }

      density.histroy = rep(NA, nchain)
      density.histroy[1] = log.density(parm=theta.init)

      alpha.candidates = matrix(NA, nrow=nchain, ncol=J)

      for (j in 1:J){
            alpha.candidates[,j] = rgamma(n=nchain, , shape=jump.shape, scale=jump.scale)
      }
      if (J>1) alpha.candidates = rowCumsums(alpha.candidates)

      for (i in 2: nchain){
            #print(i)
            #theta.new = Jump2Next(post.chain[i-1,], J=J, iter=i, alphalist=alpha.candidates, jump.sd=jump.sd2)
            theta.new = Jump2Next.b(post.chain[i-1,], J=J, jump.sd1 = jump.sd1, jump.sd2 = jump.sd2)
            density.new = log.density(parm=theta.new)
            MH.ratio = exp( density.new - density.histroy[i-1] )

            if (MH.ratio >1){
                  # accept
                  post.chain[i, ] = theta.new
                  density.histroy[i] = density.new
                  acceptance[i] = 1
            } else {
                  r = runif(1,0,1)
                  if (r < MH.ratio) {
                        # accept
                        post.chain[i, ] = theta.new
                        density.histroy[i] = density.new
                        acceptance[i] = 1
                  } else {
                        # reject
                        post.chain[i, ] = post.chain[i-1,]
                        density.histroy[i] = density.histroy[i-1]
                  }
            }
            if (i %% 26==0) {
                  print(c(acceptance[i], post.chain[i, ]))
            }
            if (i %% 50==0) {
                  print(paste("iter: ",i, ", avg. acceptance: ", round(mean(acceptance[2:i]) , 3), 
                        " posterior mean: " ))
                  print(apply(post.chain[2:i,], 2, mean))
            }
      }
      return( list(theta.posterior=post.chain, acceptance=acceptance, log.density=density.histroy) )
}




















