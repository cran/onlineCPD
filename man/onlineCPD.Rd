\name{onlineCPD}
\alias{onlineCPD}
\title{
Online Changepoint Detection
}
\description{
Performs the oCPD algorithm in a truely on-line fashion; each iteration, you provide it with one more data point in the series, and the function will update the run length probabilities.
}
\usage{
onlineCPD(oCPD = NULL, datapt, timept = NULL, hazard_func = const_hazard)
}
\arguments{
  \item{oCPD}{
an object of type "oCPD", \code{datapt} will be added to the end of this series
}
  \item{datapt}{
the next data point in the series. This must have the same number of variables as \code{oCPD}
}
  \item{timept}{
optional, the corresponding time for \code{datapt}
}
  \item{hazard_func}{
hazard function used in the model. Defaults to a constant hazard, suitable for exponential family models.
}
}
\details{
This function performs the oCPD algorithm, one point at a time. To use, start with either the result to a call to \code{\link{offlineCPD}} or just a single data point. In the case of a single data point, leave argument \code{oCPD} as \code{NULL}. See Examples section. 

See \code{\link{offlineCPD}} for more information on how the algorithm works.
}
\value{
An object of class "oCPD", which is a list containing the following:
\item{R}{n by n matrix of run-length probabilities. The value at \code{R[i,j]} is the probability that at data point \code{j}, the current run length is \code{i}.}
\item{data}{same as the input parameter, included for plotting.}
\item{time}{same as the input parameter, included for plotting.}
\item{alpha}{the vector of values of alpha after the final data point.}
\item{beta}{the vector (or matrix) of values of beta (the variance) after the final data point.}
\item{kappa}{the vector of values of kappa after the final data point.}
\item{mu}{the vector (or matrix) of values of mu (the mean) after the final data point}
\item{max}{vector of values; \code{max[i]} is the runlength with the highest probability. Used to plot the red diagonal line in \link{plot.oCPD}}
\item{changes}{locations of detected changepoints. When the algorithm can not detect the exact location of a change, multiple possible values are reported.}
}
\references{
Adams, R. P. and Mackay, D. J. C. (2007), \href{http://arxiv.org/pdf/0710.3742.pdf}{Bayesian Online Changepoint Detection}

####OUR PAPER, WHEN IT EXISTS
}
\author{
Zachary Zanussi
}
\note{
For a version of this same algorithm that functions offline (i.e. you give it the whole time series at once), see \code{\link{offlineCPD}}

This version runs a little slower because it needs to copy the "oCPD" object with every iteration. Use \code{\link{offlineCPD}} unless you don't have the whole time series available at the start.
}
\examples{
##### Univariate Data #####
set.seed(6)
x <- c(rnorm(50,mean=0.3,sd=0.15),rnorm(40,mean=0.7,sd=0.1),rnorm(60,mean=0.5,sd=0.15))
res <- onlineCPD(datapt=x[1])
for(k in x)
  res <- onlineCPD(res,k)
plot(res)

##### Real Multivariate Data #####
data(WalBelSentiment)
data(WalBelTimes)
res <- onlineCPD(datapt=WalBelSentiment[1400,],timept=WalBelTimes[1400])
for(k in 1401:1600)
  res <- onlineCPD(res,WalBelSentiment[k,],WalBelTimes[k])
plot(res)

## You can use onlineCPD to add points to an existing "oCPD" object
y <- c(rnorm(50,0.5,0.1),rnorm(20,0.48,0.02),rnorm(50,0.5,0.1))
res <- offlineCPD(y)
plot(res)
x <- rnorm(75,0.7,0.4)
for(k in x)
  res <- onlineCPD(res,k)
plot(res)
}
\keyword{ ~changepoint }
\keyword{ ~online }
