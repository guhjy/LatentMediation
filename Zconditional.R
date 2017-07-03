#### Z conditional distribution
  Zconditional <- function(y, m, x, w1, w2, 
                           UU1,UU2,UU3,UU4,UU5,UU6,
                           theta,piZ,piU1,piU2,piU3,piU4,piU5,piU6,
                           sigma1,sigma2,KNOT){
    
    piK_pre <- NULL
    U3 <- rep(0,dim(piU3)[2]); U3[UU3] <- 1
    U4 <- rep(0,dim(piU4)[2]); U4[UU4] <- 1
    U5 <- rep(0,dim(piU5)[2]); U5[UU5] <- 1
    U6 <- rep(0,dim(piU6)[2]); U6[UU6] <- 1
    
    for(k in 1:K){
      M_K <- (abs(outer(m,KNOT[,k],"-")))^3
      OMEGA <- (abs(outer(KNOT[,k],KNOT[,k],"-")))^3
      svd.OMEGA <- svd(OMEGA)
      sqrt.OMEGA <- t(svd.OMEGA$v %*% (t(svd.OMEGA$u)*sqrt(svd.OMEGA$d)))
      Mknot <- t(solve(sqrt.OMEGA, t(M_K)))
      
      piK_pre[k] <- dnorm(m, theta[1:4,k]%*%c(1,x,w1,w2), sqrt(sigma1[k]))*
        dnorm(y,theta[5:(9+knot),k]%*%c(1,x,m,Mknot,w1,w2), sqrt(sigma2[k]))*
        piZ[k]*dnorm(UU1, piU1[k,1], 1/sqrt(piU1[k,2]))*
        dnorm(UU2, piU2[k,1], 1/sqrt(piU2[k,2]))*
        dmultinom(U3, prob=piU3[k,])*
        dmultinom(U4, prob=piU4[k,])*
        dmultinom(U5, prob=piU5[k,])*
        dmultinom(U6, prob=piU6[k,])
    }
    piK <- piK_pre/sum(piK_pre)
    piK <- ifelse(is.na(piK), rep(1/K,K), piK)
    return(piK)
  }
  