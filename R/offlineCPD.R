offlineCPD <-
function(data, time = NULL, hazard_func = const_hazard,
                      m=0, k=0.01, a=0.01, b=0.0001){
  if(!is.matrix(data)) {
    if(is.null(dim(data)[2])) data <- matrix(data,ncol=1)
    else                      data <- as.matrix(data)
  }
  lambda <- 2000
  T      <- dim(data)[1]
  dim    <- dim(data)[2]
  
  muT <- mu0 <- matrix(rep(m,times=dim),ncol=dim); kappaT <- kappa0 <- k; 
  betaT <- beta0 <- matrix(rep(b,times=dim),ncol=dim); alphaT <- alpha0 <- a; 
  R      <- matrix(0,nrow=T+1,ncol=T+1)
  R[1,1] <- 1
  
  maxes  <- vector(mode="integer",length = T + 1)
  cps    <- vector(mode="integer",length = T + 1)
  
  for(t in 1:T){
    predProbs <- studentpdf(data[t,],muT,betaT*(kappaT+1)/(alphaT*kappaT),2*alphaT)
    
    H <- hazard_func(t,lambda)
    
    R[2:(t+1),(t+1)] <- R[1:t,t] * apply(as.matrix(predProbs),1,prod) * (1 - H)
    R[1,t+1]         <- sum( R[1:t,t] * apply(as.matrix(predProbs),1,prod) * H )
    
    R[,t+1]   <- R[,t+1] / sum(R[,t+1])
    
    tempmu    <- rbind(mu0,t(t(kappaT*muT) + data[t,]) / (kappaT+1))
    tempkappa <- append(kappa0,(kappaT + 1))
    tempalpha <- append(alpha0,(alphaT + 0.5))
    tempbeta  <- rbind(beta0,(betaT + (kappaT * t(data[t,]-t(muT))^2)/(2*(kappaT+1))))
    muT       <- tempmu
    kappaT    <- tempkappa
    alphaT    <- tempalpha
    betaT     <- tempbeta
    
    maxes[t]  <- match(max(R[,t]),R[,t])
    cps[t]    <- t - maxes[t]
  }
  cps <- sort(unique(cps))

  
  result <- list(R=R,data=data,time=time,alpha=alphaT,beta=betaT,kappa=kappaT,mu=muT,max=maxes,changes=cps)
  class(result) <- "oCPD"
  
  return(result)
  ####SHOULD WE TRUNCATE SMALL VALUES IN R TO SAVE SPACE???
}
