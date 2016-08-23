\name{findCP}
\alias{findCP}
\title{
Wisely Pick Changepoints
}
\description{
Function that helps narrow down possible changepoints. Sometimes the \code{\link{onlineCPD}} and \code{\link{offlineCPD}} functions report an excessive number of changepoints. This function is called inside \code{\link{plot.oCPD}} and removes them based on criteria described in Details.
}
\usage{
findCP(oCPD, buffer = 10)
}
\arguments{
  \item{oCPD}{
results of a call to \code{\link{onlineCPD}} or \code{\link{offlineCPD}}, this is the time series that the function will choose changepoints for.
}
  \item{buffer}{
a whole number. Any time there are two changepoints within \code{buffer} data points, \code{findCP} will choose one to replace them both.
}
}
\details{
The function chooses changepoints based on the following criteria; first, if there are two changepoints within \code{buffer} data points, it will pick one to replace them both, based on the sum of the probabilities for each change. Then, if the algorithm ever changes it's mind on a changepoint (i.e. the red, diagonal line ever jumps up) the function will remove the changepoint.
}
\value{
A vector of changepoints.
}
\author{
Zachary Zanussi
}



\seealso{
\code{\link{plot.oCPD}}
}
\examples{
## Two changepoints within \code{buffer} points:
set.seed(5)
y <- c(rnorm(15,0.5,0.01),0.505,0.51,0.015,0.52,rnorm(15,0.6,0.05))
res <- offlineCPD(y)
plot(res,buffer=10)     ## Will replace the two changepoints
plot(res,cleanCP=FALSE) ## Will not replace the changepoints

## A changepoint is forgotten:
set.seed(5)
y <- c(rnorm(50,0.5,0.1),rnorm(20,0.48,0.02),rnorm(50,0.5,0.1),rnorm(75,0.7,0.4))
res <- offlineCPD(y)
plot(res)               ## Will remove the changepoint that has been second-guessed
plot(res,cleanCP=FALSE) ## Will not remove the changepoint
}
