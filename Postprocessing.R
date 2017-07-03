#### required libraries  
  library(lpSolve)
  library(mcclust)
  library(cluster)

#### construct the posterior similarity matrix        
  mat <-para[[500]]$Z  #<-- burn-in
  for(ii in 1:250){
      i <- ii*2+500  #<-- thinning
      mat <- rbind(mat,para[[i]]$Z)
  }  
  sim<-comp.psm(mat)
    
  
#### find the optimal partition  
  best.sim<-NULL  
  for(i in 1:9){
    best.sim[i] <- pam(1-sim, i, diss=TRUE)$silinfo$avg.width
  }  
  num.clus <- which.max(best.sim)
  best <- pam(1-sim, num.clus, diss=TRUE)
  

  
#### calculate cluster-specific effects (for the 4 clusters in the simulation study)
  
  DE <- IE <- TE <- matrix(ncol=4, nrow=250)
  pC <- matrix(ncol=4, nrow=250)
  
  # determine clusters
  x<-best$clustering[1:50]
  index11<-as(names(which.max(table(x))),mode(x))  
  x<-best$clustering[51:150]
  index22<-as(names(which.max(table(x))),mode(x))  
  x<-best$clustering[151:250]
  index33<-as(names(which.max(table(x))),mode(x))  
  x<-best$clustering[251:300]
  index44<-as(names(which.max(table(x))),mode(x))
  
  # indices for 4 clusters
  index1<-which(best$clustering==index11)
  index2<-which(best$clustering==index22)
  index3<-which(best$clustering==index33)
  index4<-which(best$clustering==index44)
  
  
  # esimate the effects
  for(ii in 1:250){
    
    i <- ii*2+500
    Z <- para[[i]]$Z
        
    KNOT.set <- para[[t]]$KNOT
    
    THETA <- para[[t]]$theta
    
    # cluster-specific direct effects
    DE[ii,1] <- mean(sapply(index1, function(x) THETA[6,Z[x]]))
    DE[ii,2] <- mean(sapply(index2, function(x) THETA[6,Z[x]]))
    DE[ii,3] <- mean(sapply(index3, function(x) THETA[6,Z[x]]))
    DE[ii,4] <- mean(sapply(index4, function(x) THETA[6,Z[x]]))
    
    
    # cluster-specific indirect effects
    IE1.temp <- NULL
    aa <-0    
    for(tt in index1){      
      MMM1 <- c(cbind(1,1,x1,x2)[index1,]%*%THETA[1:4,Z[tt]])
      MMM0 <- c(cbind(1,0,x1,x2)[index1,]%*%THETA[1:4,Z[tt]])      
      KNOT <- KNOT.set[,Z[tt]]      
      OMEGA <- (abs(outer(KNOT,KNOT,"-")))^3
      svd.OMEGA <- svd(OMEGA)
      sqrt.OMEGA <- t(svd.OMEGA$v %*% (t(svd.OMEGA$u)*sqrt(svd.OMEGA$d)))      
      MMM1_K <- (abs(outer(MMM1,KNOT,"-")))^3
      MMM1knot <- t(solve(sqrt.OMEGA, t(MMM1_K)))      
      MMM0_K <- (abs(outer(MMM0,KNOT,"-")))^3
      MMM0knot <- t(solve(sqrt.OMEGA, t(MMM0_K)))
      aa <- aa + 1
      IE1.temp[aa] <- mean(THETA[7,Z[tt]]*THETA[2,Z[tt]]+MMM1knot%*%THETA[8:12,Z[tt]]-MMM0knot%*%THETA[8:12,Z[tt]])
    }
    IE[ii,1] <- mean(IE1.temp)
        
    IE2.temp <- NULL
    aa <-0
    for(tt in index2){
      MMM1 <- c(cbind(1,1,x1,x2)[index2,]%*%THETA[1:4,Z[tt]])
      MMM0 <- c(cbind(1,0,x1,x2)[index2,]%*%THETA[1:4,Z[tt]])
      KNOT <- KNOT.set[,Z[tt]]
      OMEGA <- (abs(outer(KNOT,KNOT,"-")))^3
      svd.OMEGA <- svd(OMEGA)
      sqrt.OMEGA <- t(svd.OMEGA$v %*% (t(svd.OMEGA$u)*sqrt(svd.OMEGA$d)))
      MMM1_K <- (abs(outer(MMM1,KNOT,"-")))^3
      MMM1knot <- t(solve(sqrt.OMEGA, t(MMM1_K)))
      MMM0_K <- (abs(outer(MMM0,KNOT,"-")))^3
      MMM0knot <- t(solve(sqrt.OMEGA, t(MMM0_K)))
      aa <- aa + 1
      IE2.temp[aa] <- mean(THETA[7,Z[tt]]*THETA[2,Z[tt]]+MMM1knot%*%THETA[8:12,Z[tt]]-MMM0knot%*%THETA[8:12,Z[tt]])
    }
    IE[ii,2] <- mean(IE2.temp)
    
    
    
    IE3.temp <- NULL
    aa <-0
    for(tt in index3){  
      MMM1 <- c(cbind(1,1,x1,x2)[index3,]%*%THETA[1:4,Z[tt]])
      MMM0 <- c(cbind(1,0,x1,x2)[index3,]%*%THETA[1:4,Z[tt]])
      KNOT <- KNOT.set[,Z[tt]]
      OMEGA <- (abs(outer(KNOT,KNOT,"-")))^3
      svd.OMEGA <- svd(OMEGA)
      sqrt.OMEGA <- t(svd.OMEGA$v %*% (t(svd.OMEGA$u)*sqrt(svd.OMEGA$d)))
      MMM1_K <- (abs(outer(MMM1,KNOT,"-")))^3
      MMM1knot <- t(solve(sqrt.OMEGA, t(MMM1_K)))
      MMM0_K <- (abs(outer(MMM0,KNOT,"-")))^3
      MMM0knot <- t(solve(sqrt.OMEGA, t(MMM0_K)))
      aa <- aa + 1
      IE3.temp[aa] <- mean(THETA[7,Z[tt]]*THETA[2,Z[tt]]+MMM1knot%*%THETA[8:12,Z[tt]]-MMM0knot%*%THETA[8:12,Z[tt]])
    }
    IE[ii,3] <- mean(IE3.temp)
    
    IE4.temp <- NULL
    aa <-0
    for(tt in index4){  
      MMM1 <- c(cbind(1,1,x1,x2)[index4,]%*%THETA[1:4,Z[tt]])
      MMM0 <- c(cbind(1,0,x1,x2)[index4,]%*%THETA[1:4,Z[tt]])
      KNOT <- KNOT.set[,Z[tt]]
      OMEGA <- (abs(outer(KNOT,KNOT,"-")))^3
      svd.OMEGA <- svd(OMEGA)
      sqrt.OMEGA <- t(svd.OMEGA$v %*% (t(svd.OMEGA$u)*sqrt(svd.OMEGA$d)))
      MMM1_K <- (abs(outer(MMM1,KNOT,"-")))^3
      MMM1knot <- t(solve(sqrt.OMEGA, t(MMM1_K)))
      MMM0_K <- (abs(outer(MMM0,KNOT,"-")))^3
      MMM0knot <- t(solve(sqrt.OMEGA, t(MMM0_K)))
      aa <- aa + 1
      IE4.temp[aa] <- mean(THETA[7,Z[tt]]*THETA[2,Z[tt]]+MMM1knot%*%THETA[8:12,Z[tt]]-MMM0knot%*%THETA[8:12,Z[tt]])
    }
    IE[ii,4] <- mean(IE4.temp)
    
    # cluster-specific total effects
    TE[ii,1] <- DE[ii,1]+IE[ii,1]
    TE[ii,2] <- DE[ii,2]+IE[ii,2]
    TE[ii,3] <- DE[ii,3]+IE[ii,3]
    TE[ii,4] <- DE[ii,4]+IE[ii,4]
    
    
    pC[ii,1] <- length(index1)/size
    pC[ii,2] <- length(index2)/size
    pC[ii,3] <- length(index3)/size
    pC[ii,4] <- length(index4)/size
  }
  
  TE<-apply(TE, 2, mean)
  IE<-apply(IE, 2, mean)
  DE<-apply(DE, 2, mean)
  pC<-apply(pC, 2, mean)
  
  


