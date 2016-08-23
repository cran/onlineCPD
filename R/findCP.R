findCP <- function(oCPD,buffer=10) {
  if(class(oCPD)!="oCPD") stop("argument oCPD must be of type \"oCPD\"")
  if(buffer > length(oCPD$max)) stop("buffer must be less than the number of data points")
  
  max <- oCPD$max
  imax <- vector("numeric",length(max))
  for(i in 1:length(max)) imax[i] <- i - max[i]
  
  changes <- sort(unique(imax))
  while(any(diff(changes)<=buffer)){
    changes[c(k1 <- which(diff(changes)<=buffer),k2 <- k1 + 1)]
    for(i in 1:length(k1)){
      tiedpair <- c(k1[i],k2[i])
      tiedrun <- which(imax == changes[tiedpair[1]] | imax == changes[tiedpair[2]])
      v1 <- min(changes[tiedpair[1]],changes[tiedpair[2]]); v2 <- max(changes[tiedpair[1]],changes[tiedpair[2]])
      difference <- v2 - v1
      minran <- v1 - ceiling(buffer / 2) + ceiling(difference / 2); maxran <- v2 + ceiling(buffer / 2) - ceiling(difference / 2)
      tiedrun <- tiedrun[tiedrun > maxran]
      runProbs <- sapply(minran:maxran,function(val) return(sum(diag(oCPD$R[tiedrun-val,tiedrun]))))
      imax[which(imax == changes[tiedpair[1]] | imax == changes[tiedpair[2]])]  <-  minran + which.max(runProbs)
    }
    changes <- sort(unique(imax))
  }
  
  val <- imax[length(imax)-1]

  for(i in length(imax):1) {
    if(imax[i] > val) {
      imax[i] <- 0
    } else {
      val <- imax[i]
    }
  } 
  return(unique(imax))
}