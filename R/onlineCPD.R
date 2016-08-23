onlineCPD <- 
function(oCPD=NULL,datapt,timept=NULL,hazard_func=const_hazard) {
  if(is.null(oCPD)) {
    return(offlineCPD(t(datapt),timept))
  }
  if(!(class(oCPD)=="oCPD")) stop("argument oCPD must be of type \"oCPD\"")
    
  if(dim(oCPD$data)[2]!=length(datapt)) stop("oCPD and datapt must have the same number of variables")
  
  T <- dim(oCPD$R)[1]
  R2 <- matrix(0,nrow=T+1,ncol=T+1)
  R2[1:(T),1:(T)] <- oCPD$R
  lambda <- 2000
  
  predProbs <- studentpdf(datapt,oCPD$mu,oCPD$beta*(oCPD$kappa+1)/(oCPD$alpha*oCPD$kappa),2*oCPD$alpha)
  
  H <- hazard_func(T,lambda)
  
  R2[2:(T+1),(T+1)] <- R2[1:T,T] * apply(as.matrix(predProbs),1,prod) * (1 - H)
  R2[1,T+1]         <- sum( R2[1:T,T] * apply(as.matrix(predProbs),1,prod) * H )
  
  R2[,T+1]   <- R2[,T+1] / sum(R2[,T+1])
  
  tempmu    <- rbind(0,t(t(oCPD$kappa*oCPD$mu) + datapt) / (oCPD$kappa+1))
  tempkappa <- append(0.01,(oCPD$kappa + 1))
  tempalpha <- append(0.01,(oCPD$alpha + 0.5))
  tempbeta  <- rbind(0.0001,(oCPD$beta + (oCPD$kappa * t(datapt-t(oCPD$mu))^2)/(2*(oCPD$kappa+1))))
  
  maxes <- append(oCPD$max,match(max(R2[,T+1]),R2[,T+1]))
  cps   <- append(oCPD$changes,T+1 - maxes[T+1])
  cps   <- sort(unique(cps))
  
  result <- list(R=R2,data=rbind(oCPD$data,datapt,deparse.level=0),time=append(oCPD$time,timept),
                 alpha=tempalpha,beta=tempbeta,kappa=tempkappa,mu=tempmu,max=maxes,changes=cps)
  class(result) <- "oCPD"
  return(result)
}