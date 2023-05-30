## MCMC Gibbs and MH algorithm for spatial R2D2 prior

# Y: response variable
# X: explanatory variables
# Z: random effects design matrix
# s: spatial sampling locations
# nu: shape parameter for Matern model
# a,b: Prior distribution for R2 ~ Beta(a,b)
# a0, b0: Hyper-paramters for sigma2 ~ IG(a0, b0)
# mu0, sigma20: Hyper-parameter for beta0 ~ N(mu0, sigma20)
# xi0: Hyper-parameter for phi ~ Dirichlet(xi0, ..., xi0)
# MH_spat: Standard deviation for MH proposal values for rho
# cphi: Multiplier for MH proposal values for phi
# iters: number of MCMC iterations
# burn: number of MCMC burn-in samples

mpa_mcmc <- function(Y, X, Z, s, nu=0.5, a=1, b=1, a0=0.1, b0=0.1, mu0=0, sigma20=100,xi0=1, MH_spat=0.2,cphi=20, iters=1000,burn=1000){
  library(mvnfast)
  library(geoR)
  library(extraDistr)
  library(GeneralizedHyperbolic)
  tick = proc.time()[3]
  
  
  # Bookkeeping
  n        <- length(Y)
  p        <- ncol(X)
  L        <- ncol(Z)
  P        <- diag(n) - matrix(1/n, nrow=n, ncol=n)
  
  tXX      <- t(X)%*%X
  tZZ      <- t(Z)%*%Z
  
  # Initial values
  beta0    <- mean(Y)
  beta     <- lm(Y ~ X - 1)$coef # fixed effects
  XB       <- X%*%beta
  U        <- rnorm(L)                 # non-spatial random effects
  ZU       <- Z%*%U
  theta    <- Y - XB - Z%*%U           # spatial random effects
  sigma2   <- as.numeric(var(Y))
  
  gmm      <- rgamma(1, b, 1)
  u        <- rgamma(1, a, gmm)
  v        <- rgamma(1, 0.1, 0.1)
  w        <- u*v                    # global variance
  
  phi      <- rep(1/3, 3)  # apportions variance between different effects
  
  logrho   <- log(0.5)
  lognu    <- log(nu)
  param    <- c(logrho, lognu)
  d        <- as.matrix(rdist(s))
  maxd     <- max(d)
  Sigma    <- matern_cov(d, exp(logrho), exp(lognu))
  SigmaInv <- solve(Sigma)
  TST      <- t(theta)%*%SigmaInv%*%theta
  
  
  PXX   <- P%*%X%*%t(X)
  PZZ   <- P%*%Z%*%t(Z)
  PSig  <- P%*%Sigma
  
  eigs    <- eigen(phi[1]/p*PXX + phi[2]*PZZ + phi[3]*PSig, TRUE, TRUE)$values
  mux     <- sum(eigs)/(n-1)
  sigma2x <- 2*sum(eigs^2)/(n-1)^2
  
  aa <- mux^2/sigma2x; bb=sigma2x/mux
  
  # Keep track of stuff
  
  keep.beta  <- matrix(0, nrow=iters, ncol=p)
  keep.theta <- matrix(0, nrow=iters, ncol=n)
  keep.b0    <- matrix(0, nrow=iters, ncol=1)
  keep.phi   <- matrix(0, nrow=iters, ncol=3)
  keepers    <- matrix(0, nrow=iters, ncol=4)
  keep.r2    <- matrix(0, nrow=iters, ncol=1)
  acc1       <- 0
  acc2       <- 0
  
  
  colnames(keepers)   <- c("sigma2","W", "rho", "nu")
  
  # GO!!!
  
  for(iter in 1:(iters+burn)){
    
    ##############################################:
    #####       MEAN PARAMETERS (Gibbs)    #######:
    ##############################################:
    
    VVV      <- (n/sigma2 + 1/sigma20)^(-1)
    MMM      <- sum(Y-XB-ZU-theta)/sigma2 + mu0/sigma20
    beta0    <- rnorm(1, MMM*VVV, sqrt(VVV))
    
    VVV      <- solve(tXX + diag(p)/(phi[1]*w/p) )
    MMM      <- t(X)%*%(Y-beta0-ZU-theta)
    beta     <- VVV%*%MMM + sqrt(sigma2)*t(chol(VVV))%*%rnorm(p)
    
    XB       <- X%*%beta
    
    VVV      <- solve(tZZ + diag(L)/(phi[2]*w))
    MMM      <- t(Z)%*%(Y-beta0-XB-theta)
    U        <- VVV%*%MMM + sqrt(sigma2)*t(chol(VVV))%*%rnorm(L)
    
    ZU       <- Z%*%U
    
    VVV      <- solve(diag(n) + SigmaInv/(phi[3]*w) )
    MMM      <- Y-beta0-XB-ZU
    theta    <- VVV%*%MMM + sqrt(sigma2)*t(chol(VVV))%*%rnorm(n)
    
    TST      <- t(theta)%*%SigmaInv%*%theta
    
    ##############################################:
    #####          VARIANCES (Gibbs)        #######:
    ##############################################:
    
    # Response
    sigma2    <- extraDistr::rinvgamma(1, alpha=a0 + n + p/2 + L/2,
                                       beta=b0 + t(Y-beta0-XB-ZU-theta)%*%(Y-beta0-XB-ZU-theta)/2 +
                                         sum(beta^2)/(2*phi[1]*w/p)+
                                         sum(U^2)   /(2*phi[2]*w)+
                                         TST        /(2*phi[3]*w))
    
    # Global
    u <- rgig(1, lambda = a - (n+p+L)/2,
              psi = 2*gmm,
              chi = sum(beta^2)/(sigma2*v*phi[1]/p) + sum(U^2)/(sigma2*v*phi[2]) +
                TST/(phi[3]*v*sigma2))
    
    v <- extraDistr::rinvgamma(1, alpha = aa + (n+p+L)/2,
                               beta = 1/bb + 0.5*sum(beta^2)/(sigma2*u*phi[1]/p) + 
                                 0.5*sum(U^2)/(sigma2*u*phi[2]) +
                                 0.5*TST/(phi[3]*u*sigma2))
    w <- u*v
    
    ##############################################:
    #####        EXTRA PARAMS (Gibbs)      #######:
    ##############################################:
    
    gmm <- rgamma(1, shape=a+b, rate=1+u)
    
    
    ##############################################:
    ######## VARIANCE WEIGHTS (Metropolis) #######:
    ##############################################:
    
    ############## PHI
    curll = sum(dnorm(beta, 0, sqrt(sigma2*phi[1]/p*w), log=T)) +
      sum(dnorm(U, 0, sqrt(sigma2*phi[2]*w), log=T))+
      dmvn(as.numeric(theta), rep(0,n), sigma2*phi[3]*w*Sigma, log=T)+
      dbetapr(w, a, aa, 1/(bb*gmm), log=T) +
      ddirichlet(phi, rep(xi0, 3), log=T)
    
    canphi <- rdirichlet(1, cphi * phi)
    
    while(max(canphi) > 0.99 | min(canphi) < 0.01){
      canphi <- rdirichlet(1, cphi * phi)
      canPHI <- diag(canphi[1:p])
    }
  
    caneigs    <- eigen(canphi[1]/p*PXX + canphi[2]*PZZ + canphi[3]*PSig, TRUE, TRUE)$values
    canmux     <- sum(caneigs)/(n-1)
    cansigma2x <- 2*sum(caneigs^2)/(n-1)^2
    canaa <- canmux^2/cansigma2x; canbb=cansigma2x/canmux
    
    canll = sum(dnorm(beta, 0, sqrt(sigma2*canphi[1]/p*w), log=T)) +
      sum(dnorm(U, 0, sqrt(sigma2*canphi[2]*w), log=T))+
      dmvn(as.numeric(theta), rep(0,n), sigma2*canphi[3]*w*Sigma, log=T)+
      dbetapr(w, a, canaa, 1/(canbb*gmm), log=T) +
      ddirichlet(canphi, rep(xi0, 3), log=T)
    
    MH <- canll -  curll +
      ddirichlet(phi, canphi, log=T) - ddirichlet(canphi, phi, log=T)
    
    
    if(log(runif(1)) < MH){
      phi     <- canphi
      aa      <- canaa
      bb      <- canbb
      
      acc1 =  acc1 + 1
    }
    
    if(iter <= burn & iter%%100==0){
      
      ar1 = acc1 / 100
      
      if( ar1 > 0.5){
        cphi <- cphi / 2
      }else if(ar1 < 0.2){
        cphi <- cphi * 2
      }
      
      acc1 <- 0
      
    }
    
    
    ##############################################:
    #### CORRELATION PARAMETERS (Metropolis) #####:
    ##############################################:
    
    # rho and nu (spatial range and smoothness parameter)

    curll <- dmvn(as.numeric(theta), rep(0, n), sigma2*phi[3]*w*Sigma, log=T ) +
      dnorm(logrho, -2, sqrt(1), log=T) +
      dbetapr(w, a, aa, 1/(bb*gmm), log=T)

    # Propose candidate value
    canparam   <- c(param[1] + MH_spat*rnorm(1), lognu)
    canSig    <- matern_cov(d, exp(canparam[1]), exp(canparam[2]) )
    canSig[is.na(canSig)] <- 1
    canSig[canSig > 1] <- 1

    PcanSig    <- P%*%canSig
    caneigs    <- eigen(phi[1]/p*PXX + phi[2]*PZZ + phi[3]*PcanSig, TRUE, TRUE)$values
    canmux     <- sum(caneigs)/(n-1)
    cansigma2x <- 2*sum(caneigs^2)/(n-1)^2
    canaa <- canmux^2/cansigma2x; canbb=cansigma2x/canmux

    canll <- dmvn(as.numeric(theta), rep(0, n), sigma2*phi[3]*w*canSig, log=T ) +
      dnorm(canparam[1], -2, sqrt(1), log=T) +
      dbetapr(w, a, canaa, 1/(canbb*gmm), log=T)

    MH  <- canll - curll

    if(log(runif(1)) < MH){
      logrho   <- canparam[1]
      param    <- c(logrho, lognu)
      Sigma    <- canSig
      SigmaInv <- solve(Sigma)
      aa       <- canaa
      bb       <- canbb
      PSig     <- PcanSig

      acc2       <- acc2 + 1
    }

    if(iter <= burn & iter%%100==0){

      ar2 = acc2 / 100

      if( ar2 > 0.5){
        MH_spat <- MH_spat * 1.2
      }else if(ar2 < 0.2){
        MH_spat <- MH_spat * 0.8
      }

      acc2 <- 0

    }
    
    if(iter > burn & iter%%1000==0){
      print(iter-burn)
    }
    
    
    ##############################################:
    #####        KEEP TRACK OF STUFF       #######:
    ##############################################:
    if(iter>burn){
      keep.b0[iter-burn,]    <- beta0
      keep.beta[iter-burn,]  <- beta
      keep.theta[iter-burn,]  <- theta
      keep.phi[iter-burn,]   <- phi
      keepers[iter-burn,]    <- c(sigma2, w, exp(logrho), exp(lognu))
      
      r2 <- var(beta0 + X%*%beta + Z%*%U + theta)
      keep.r2[iter-burn,] <- r2/(r2 + sigma2)
      
      
      if(length(pred) > 0){
        pred.hat <- beta0 + Xpr%*%beta + Zpr%*%U + dpr%*%SigmaInv%*%(Y-beta0-XB-ZU)
        keep.pred[iter-burn, ] <- t(pred.hat)
      }
      
    }
  }
  
  tock = proc.time()[3]
  

  list(beta=keep.beta, theta=keep.theta, keepers=keepers,phi=keep.phi,beta0=keep.b0,r2=keep.r2,
         minutes=(tock-tick)/60,
         ar_spat = acc2/iters,
         ar_phi  = acc1/iters,
         MH_spat=MH_spat,
         cphi = cphi)
  
}
