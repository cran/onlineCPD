\name{onlineCPD-package}
\alias{onlineCPD-package}
\alias{onlineCPDpackage}
\docType{package}
\title{
onlineCPD
}
\description{
Detects multiple changepoints in uni- or multivariate time series data. The algorithm is based on Bayesian methods and detects changes on-line; ie the model updates with every observation rather than relying on retrospective segmentation. However, the user may choose to use the algorithm off- or on-line. 
}
\details{

The DESCRIPTION file:
Detects multiple changepoints in uni- or multivariate time series data. The algorithm is based on Bayesian methods and detects changes on-line; ie the model updates with every observation rather than relying on retrospective segmentation. However, the user may choose to use the algorithm off- or on-line. 


This package provides two ways to use the oCPD algorithm; first is the function \code{\link{offlineCPD}}, where you provide your whole time series at once. The algorithm still performs as if it were being provided data points one at a time, but it's more convenient for certain analysis. Second, there is \code{\link{onlineCPD}}, where you provide your time series one point at a time, as the oCPD algorithm is meant to be performed. 

BETA Release
}
\author{
Zachary Zanussi <zachary9506@gmail.com>

Maintainer: Zachary Zanussi <zachary9506@gmail.com>
}
\references{
Adams, R. P. and Mackay, D. J. C. (2007),  \href{http://arxiv.org/pdf/0710.3742.pdf}{Bayesian Online Changepoint Detection}

###OUR PAPER,WHEN IT EXISTS
}
\keyword{ package }
\seealso{
\link{onlineCPD}
}
\examples{
##### Univariate Data #####
set.seed(6)
x <- c(rnorm(50,mean=0.3,sd=0.15),rnorm(40,mean=0.7,sd=0.1),rnorm(60,mean=0.5,sd=0.15))
plot(offlineCPD(x))

##### Real Multivariate Data #####
data(WalBelSentiment)
data(WalBelTimes)
plot(offlineCPD(WalBelSentiment[1400:1600,],WalBelTimes[1400:1600]))
}
