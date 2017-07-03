#### Main MCMC function
  MCMC <- function(piZ,gamma,
                   piU1,piU2,piU3,piU4,piU5,piU6, 
                   theta,sigma1,sigma2,chi,ppiL,sigmaf,KNOT){
    
    # upate cluster
    Z <- apply(cbind(Y,M,v,X), 1, function(x) 
      which(rmultinom(1,1,Zconditional(x[1],x[2],x[3],x[4],x[5],
                                       x[4],x[5],x[6],x[7],x[8],x[9],
                                       theta,piZ,piU1,piU2,piU3,piU4,piU5,piU6,
                                       sigma1,sigma2,KNOT))==1))
    
    # update weights
    W <- NULL
    for(c in 1:(K-1)){
      W[c] <- rbeta(1, 1+length(which(Z==c)), gamma+sum(length(which(Z > c))))
    }
    piZ <- NULL
    piZ_temp <- NULL
    piZ_temp[1] <- 1
    piZ[1] <- W[1]
    for(j in 2:(K-1)){
      piZ_temp[j] <- piZ_temp[j-1]*(1-W[j-1])
      piZ[j] <- W[j]*piZ_temp[j]
    }
    piZ[K] <- piZ_temp[K-1]*(1-W[K-1])
    
    
    # update gamma parameter
    gamma <- rgamma(1, 1+K-1, 1-sum(log(1-W)))
        
    
    # spline for the mediator term
    KNOT <- matrix(nrow=knot, ncol=K)
    for(k in 1:K){
      if(length(which(Z==k)) > knot){
        KNOT[,k] <- quantile(M[which(Z==k)], seq(0,1,length=knot+2)[-c(1,(knot+2))])
      }else{
        KNOT[,k] <- quantile(M, seq(0,1,length=knot+2)[-c(1,(knot+2))])
      }
    }
    
    # update piU1 ... piU6    
    piU3_pre <- matrix(nrow=K,ncol=3)
    piU4_pre <- matrix(nrow=K,ncol=4)
    piU5_pre <- matrix(nrow=K,ncol=4)
    piU6_pre <- matrix(1/5,nrow=K,ncol=5)
    
    for(k in 1:K){
      # update piU1
      X1.len <- length(X[which(Z==k),1])
      weight1 <- X1.len*piU1[k,2]/(X1.len*piU1[k,2]+1*piU1[k,2])
      if(X1.len > 0){
        mean1 <- mean(X[which(Z==k),1])
      }else{
        mean1 <- 0
      }
      piU1[k,1] <- rnorm(1, weight1*mean1+(1-weight1)*mean(X[,1]), sqrt(1/(X1.len*piU1[k,2]+1*piU1[k,2])))
      if(X1.len > 0){
        weight11 <- 0.5*sum((X[which(Z==k),1]-mean1)^2)+X1.len*1/(2*(X1.len+1))*(mean1-mean(X[,1]))^2
      }else{
        weight11<-0
      }
      piU1[k,2] <- rgamma(1,1+X1.len/2, 1+weight11)
      
      # update piU2
      X2.len <- length(X[which(Z==k),2])
      weight2 <-X2.len*piU2[k,2]/(X2.len*piU2[k,2]+1*piU2[k,2])
      if(X2.len > 0){
        mean2 <- mean(X[which(Z==k),2])
      }else{
        mean2 <- 0
      }
      piU2[k,1] <- rnorm(1, weight2*mean2+(1-weight2)*mean(X[,2]), sqrt(1/(X2.len*piU2[k,2]+1*piU2[k,2])))
      if(X2.len > 0){
        weight22 <- 0.5*sum((X[which(Z==k),2]-mean2)^2)+X2.len*1/(2*(X2.len+1))*(mean2-mean(X[,2]))^2
      }else{
        weight22 <- 0
      }
      piU2[k,2] <- rgamma(1,1+X2.len/2, 1+weight22)
            
      # update piU3, piU4, piU5, piU6
      piU3_pre[k,] <- table( factor(X[which(v==1),3], levels = 1:3))+1
      piU4_pre[k,] <- table( factor(X[which(v==1),4], levels = 1:4))+1
      piU5_pre[k,] <- table( factor(X[which(v==1),5], levels = 1:4))+1
      piU6_pre[k,] <- table( factor(X[which(v==1),6], levels = 1:5))+1
      
      piU3[k,] <- piU3_pre[k,]/sum(piU3_pre[k,])
      piU4[k,] <- piU4_pre[k,]/sum(piU4_pre[k,])
      piU5[k,] <- piU5_pre[k,]/sum(piU5_pre[k,])
      piU6[k,] <- piU6_pre[k,]/sum(piU6_pre[k,])
    }
    
    
    # update coefficients for the mediator model
    alpha1 <- 15
    beta1 <- 7
    B0 <- diag(c(var(M), var(M),var(M),var(M)))
    MU <- c(0,0,0,0)
    XX <- cbind(1, v, X[,1:2])
    
    for(k in 1:K){      
      if(length(which(Z==k))==0){
        theta[1:4,k] <- rmnorm(1, MU, sigma1[k]*B0)
        sigma1[k] <- 1/rgamma(1,alpha1, beta1)
      }else{
        if(length(which(Z==k))==1){
          mu <- solve(t(t(XX[which(Z==k),]))%*%t(XX[which(Z==k),])+solve(B0))%*%(t(t(XX[which(Z==k),]))%*%M[which(Z==k)]+solve(B0)%*%MU)
          theta[1:4,k] <- rmnorm(1, solve(t(t(XX[which(Z==k),]))%*%t(XX[which(Z==k),])+solve(B0))%*%(t(t(XX[which(Z==k),]))%*%M[which(Z==k)]+solve(B0)%*%MU), sigma1[k]*solve(t(t(XX[which(Z==k),]))%*%XX[which(Z==k),]+solve(B0)) )
          sigma1[k] <- 1/rgamma(1, alpha1+(length(which(Z==k)))/2, beta1+(t(t(M[which(Z==k)]))%*%t(M[which(Z==k)]) +t(MU)%*%solve(B0)%*%MU - t(mu)%*%(t(t(XX[which(Z==k),]))%*%XX[which(Z==k),]+solve(B0))%*%mu )/2 )
        }else{
          mu <- solve(t(XX[which(Z==k),])%*%XX[which(Z==k),]+solve(B0))%*%(t(XX[which(Z==k),])%*%M[which(Z==k)]+solve(B0)%*%MU)
          theta[1:4,k] <- rmnorm(1, solve(t(XX[which(Z==k),])%*%XX[which(Z==k),]+solve(B0))%*%(t(XX[which(Z==k),])%*%M[which(Z==k)]+solve(B0)%*%MU), sigma1[k]*solve(t(XX[which(Z==k),])%*%XX[which(Z==k),]+solve(B0)) )
          sigma1[k] <- 1/rgamma(1, alpha1+(length(which(Z==k)))/2, beta1+(t(M[which(Z==k)])%*%M[which(Z==k)] +t(MU)%*%solve(B0)%*%MU - t(mu)%*%(t(XX[which(Z==k),])%*%XX[which(Z==k),]+solve(B0))%*%mu )/2 )
        }
      }
    }
    
    # update coefficients for the outcome model    
    alpha2<-15
    beta2<-7
    
    for(k in 1:K){
      if(length(which(Z==k))==0){
        mu <- rep(0,knot+5)
        B1 <- diag(var(Y),5+knot)
        theta[5:(9+knot),k] <- rmnorm(1, mu, sigma2[k]*B1)
        sigma2[k] <- 1/rgamma(1,alpha2, beta2)
      }else{
        mu <- c(0,chi[k],0, rep(0,knot+2))
        B1 <- diag(c(var(Y), s^2/sigma2[k],var(Y),rep(sigmaf[k],knot),var(Y),var(Y)))
        OMEGA <- (abs(outer(KNOT[,k],KNOT[,k],"-")))^3
        svd.OMEGA <- svd(OMEGA)
        sqrt.OMEGA <- t(svd.OMEGA$v %*% (t(svd.OMEGA$u)*sqrt(svd.OMEGA$d)))
        
        M_K <- (abs(outer(M[which(Z==k)],KNOT[,k],"-")))^3
        Mknot <- t(solve(sqrt.OMEGA, t(M_K)))
        
        
        if(length(which(Z==k))==1){
          XX<-matrix(c(1, v[which(Z==k)], M[which(Z==k)], Mknot, X[which(Z==k),1:2]),nrow=1)
          mu_n <- solve(t(XX)%*%XX+solve(B1))%*%((t(XX))%*%Y[which(Z==k)]+solve(B1)%*%mu)
          theta[5:(9+knot),k] <- rmnorm(1, solve((t(XX))%*%XX+solve(B1))%*%((t(XX))%*%Y[which(Z==k)]+solve(B1)%*%mu), sigma2[k]*solve((t(XX))%*%XX+solve(B1)) )
          sigma2[k] <- 1/rgamma(1, alpha2+length(which(Z==k))/2, beta2+ (t(t(Y[which(Z==k)]))%*%t(Y[which(Z==k)])+t(mu)%*%solve(B1)%*%mu-t(mu_n)%*%((t(XX))%*%XX+solve(B1))%*%mu_n)/2 )
        }else{
          XX<-cbind(1, v[which(Z==k)], M[which(Z==k)], Mknot, X[which(Z==k),1:2])
          mu_n <- solve(t(XX)%*%XX+solve(B1))%*%(t(XX)%*%Y[which(Z==k)]+solve(B1)%*%mu)
          theta[5:(9+knot),k] <- rmnorm(1, solve(t(XX)%*%XX+solve(B1))%*%(t(XX)%*%Y[which(Z==k)]+solve(B1)%*%mu), sigma2[k]*solve(t(XX)%*%XX+solve(B1)) )
          sigma2[k] <- 1/rgamma(1, alpha2+length(which(Z==k))/2, beta2+ (t(Y[which(Z==k)])%*%Y[which(Z==k)]+t(mu)%*%solve(B1)%*%mu-t(mu_n)%*%(t(XX)%*%XX+solve(B1))%*%mu_n)/2)
        }
      }
    }
        
    # update the parameter of the informative prior
    piL_pre <- matrix(NA,K,K)
    piL <- NULL
    for(k in 1:K){
      for(l in 1:K){
        piL_pre[k,l] <- dnorm(theta[6,k], List.Values[l], s)*ppiL[k,l]
      }
      piL <- piL_pre[k,]/sum(piL_pre[k,])
      piL <- ifelse(is.na(piL), rep(1/K,K), piL)
      indicator <- rmultinom(1,1,piL)
      chi[k] <- List.Values[which(indicator==1)]
      ppiL[k,] <- ppiL[k,]+indicator
    }
    
    # update sigma_f 
    prop<-NULL
    for(l in 1:K){
      prop[l] <- runif(1, max(sigmaf[l]-1, 0.001), sigmaf[l]+1) 
      rat <- dmnorm(theta[8:(7+knot),l],rep(0,knot),diag(sigma2[l]*prop[l],knot,knot),log=TRUE)+
        dgamma(prop[l],20,0.5,log=TRUE)+dunif(sigmaf[l], max(prop[l]-1, 0.001), prop[l]+1, log=TRUE)-
        dmnorm(theta[8:(7+knot),l],rep(0,knot),diag(sigma2[l]*sigmaf[l],knot,knot),log=TRUE)-
        dgamma(sigmaf[l],20,0.5,log=TRUE)-dunif(prop[l], max(sigmaf[l]-1, 0.001), sigmaf[l]+1, log=TRUE)
      if (is.na(rat)){
        prop[l] <- sigmaf[l]
      }else{
        if (log(runif(1))>rat) {
          prop[l] <- sigmaf[l]
        }
      }
    }
    
    return(list(piZ=piZ,
                gamma=gamma,
                piU1=piU1,
                piU2=piU2,
                piU3=piU3,
                piU4=piU4,
                piU5=piU5,
                piU6=piU6, 
                theta=theta,
                sigma1=sigma1,
                sigma2=sigma2,
                chi=chi,
                ppiL=ppiL,
                KNOT=KNOT,
                prop=prop,
                Z=Z))
  }
 
  