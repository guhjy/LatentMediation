#### generate data
  source("SampleData.R")

#### import functions
  source("Zconditional.R")
  source("MainMCMC.R")

#### required libraries
  library(mnormt)
  
#### the numbers of clusters & knots
  K <- 21
  knot <- 5

#### set the range of lambda 
  fit <- lm(Y ~ v+X[,1:2]) #<------- check this
  List.Values <- seq(coef(fit)[2]-sqrt(var(Y[which(v==1)])+var(Y[which(v==0)])),coef(fit)[2]+sqrt(var(Y[which(v==1)])+var(Y[which(v==0)])),length=K)
  s <- (List.Values[2]-List.Values[1])/2  # distanced by 2 standard deviation.
  
  
#### initials
  piZ.init <- rep(1,K)/K
  gamma.init <- 2  
  chi2.init <- rep(0,K)
  ppiL2.init <- matrix(1, K,K)
  
  piU1.init <- cbind(rep(0,K),rep(1,K))
  piU2.init <- cbind(rep(0,K),rep(1,K))
  piU3.init <- matrix(1/3,nrow=K,ncol=3)
  piU4.init <- matrix(1/4,nrow=K,ncol=4)
  piU5.init <- matrix(1/4,nrow=K,ncol=4)
  piU6.init <- matrix(1/5,nrow=K,ncol=5)
  theta.init <- matrix(0, nrow=9+knot, ncol=K)
  sigma1.init <- rep(var(M,na.rm=TRUE),K)
  sigma2.init <- rep(var(Y,na.rm=TRUE),K)
  chi.init <- rep(0,K)
  ppiL.init <- matrix(1, K,K)
  sigmaf.init <- rep(var(Y,na.rm=TRUE),K)
  
  KNOT.init <- quantile(M, seq(0,1,length=knot+2)[-c(1,(knot+2))])
  KNOT.init <- matrix(rep(KNOT.init,K), ncol=K, byrow=FALSE)
      
  
  MCMC.num <- 1000 #<--- the number of iterations
  para <- list()
  para[[1]] <- list(piZ=piZ.init,
                    gamma=gamma.init,
                    piU1=piU1.init,
                    piU2=piU2.init,
                    piU3=piU3.init,
                    piU4=piU4.init,
                    piU5=piU5.init,
                    piU6=piU6.init, 
                    theta=theta.init,
                    sigma1=sigma1.init,
                    sigma2=sigma2.init,
                    chi=chi.init,
                    ppiL=ppiL.init,
                    KNOT=KNOT.init,
                    prop=sigmaf.init,
                    Z=rep(NA,size))
  
#### MCMC Iterations
  pb <- txtProgressBar(min = 0, max = MCMC.num, style = 3)
  for (t in 2:MCMC.num){
            
    para[[t]] <- MCMC(piZ = para[[t-1]]$piZ,
                     gamma = para[[t-1]]$gamma,
                     piU1 = para[[t-1]]$piU1,
                     piU2 = para[[t-1]]$piU2,
                     piU3 = para[[t-1]]$piU3,
                     piU4 = para[[t-1]]$piU4,
                     piU5 = para[[t-1]]$piU5,
                     piU6 = para[[t-1]]$piU6,
                     theta = para[[t-1]]$theta,
                     sigma1 = para[[t-1]]$sigma1,
                     sigma2 = para[[t-1]]$sigma2,
                     chi = para[[t-1]]$chi,
                     ppiL = para[[t-1]]$ppiL,
                     sigmaf = para[[t-1]]$prop,
                     KNOT =  para[[t-1]]$KNOT
    )
    print(t)
  }
  
  
