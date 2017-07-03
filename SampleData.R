#### required libraries
  library(sn)

#### sample size = 300
  size <- 300

#### 3 clusters by the covariates level
  clus1 <- c(rep(1,50),rep(2,200),rep(3,50))

#### generate covariates
  X <- matrix(nrow=size,ncol=6)
  
  for(k in 1:size){
    if(clus1[k]==1){
      X[k,1] <- rnorm(1,0,0.1)
      X[k,2] <- rnorm(1,1,0.1)
      X[k,3] <- which(rmultinom(1, 1, c(0.1,0.2,0.7))==1)
      X[k,4] <- which(rmultinom(1, 1, c(0,0.1,0.1,0.8))==1)
      X[k,5] <- which(rmultinom(1, 1, c(0,0.1,0.3,0.6))==1)
      X[k,6] <- which(rmultinom(1, 1, c(0,0,0.1,0.2,0.7))==1)
    }
    
    if(clus1[k]==2){
      X[k,1] <- rnorm(1,1,0.1)
      X[k,2] <- rnorm(1,2,0.1)
      X[k,3] <- which(rmultinom(1, 1, c(0.7,0.2,0.1))==1)
      X[k,4] <- which(rmultinom(1, 1, c(0.8,0.1,0.1,0))==1)
      X[k,5] <- which(rmultinom(1, 1, c(0.6,0.3,0.1,0))==1)
      X[k,6] <- which(rmultinom(1, 1, c(0.7,0.2,0.1,0,0))==1)
    }
    
    if(clus1[k]==3){
      X[k,1] <- rnorm(1,1.5,0.1)
      X[k,2] <- rnorm(1,0.5,0.1)
      X[k,3] <- which(rmultinom(1, 1, c(0.1,0.2,0.7))==1)
      X[k,4] <- which(rmultinom(1, 1, c(0.6,0.2,0.1,0.1))==1)
      X[k,5] <- which(rmultinom(1, 1, c(0.7,0.1,0.1,0.1))==1)
      X[k,6] <- which(rmultinom(1, 1, c(0,0,0.1,0.1,0.8))==1)
    }
    
    
  }
  
  
  
#### define expit() function
  expit <- function(x){
    exp(x)/(1+exp(x))
  }
  
#### Treatment Assignment Mechanism
  B <- c(0.6, 0.6, -0.5)
  Xps <- cbind(1, X[,1:2])
  pG <- expit(Xps%*%B)
  v <- rbinom(size, 1, pG)
  
  
#### 3 clusters by the outcome-mediator level  
  clus2 <- c(rep(1,50),rep(2,100),rep(3,150))

#### generate outcomes and mediators  
  Y<-NULL
  M<-NULL
  
  for(k in 1:size){
    if(clus2[k]==1){
      M[k] <- rnorm(1, 3+0.3*v[k]+0.5*X[k,1]+1*X[k,2], sd=0.5)
      Y[k] <- rnorm(1, -2-1*v[k]+0.5*sin(M[k]*10)+0.2*X[k,1]+0.4*X[k,2], sd=0.5)
    }
    
    if(clus2[k]==2){
      M[k] <- 1-0.5*v[k]+1.5*X[k,1]+1.5*X[k,2]+rnorm(1,0,0.5)
      Y[k] <- rnorm(1, -1-1*v[k]-1*M[k]+0.5*X[k,1]+0.5*X[k,2], sd=1.5)
    }
    
    if(clus2[k]==3){
      M[k] <- 3+0.6*v[k]+0.5*X[k,1]+0.5*X[k,2]+rnorm(1,0,0.5)
      Y[k] <- rsn(1, xi=1*v[k]+2*M[k]+0.3*X[k,1]+0.3*X[k,2], omega=1, alpha=3)
    }        
  }
  

  


