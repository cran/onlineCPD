\name{offlineCPD}
\alias{offlineCPD}
\title{
Offline Bayesian Changepoint Detection
}
\description{
An algorithm for detecting multiple changepoints in uni- or multivariate time series. The algorithm works on-line; ie the model is calculated and updated with each data observation. Though the algorithm performs as if data was supplied on-line, this version of the algorithm takes the whole series at once, ie it performs off-line. See \code{\link{onlineCPD}} for a version that runs iteratively, one data point at a time. The algorithm implements the Bayesian methods given in Adams and Mackay (2007) and is based on Matlab code released with the paper. The model has been extended to work on multivariate data.
}
\usage{
offlineCPD(data, time = NULL, hazard_func = const_hazard,
 m = 0, k = 0.01, a = 0.01, b = 1e-04)
}
\arguments{
  \item{data}{
a vector (for univariate) or matrix (for multivariate) composed of time series data. For multivariate, each column is a different time series. Column names will be extracted for plotting, so name the columns accordingly. Note that you must provide the whole series; if you want to provide one data point at a time, use \code{\link{onlineCPD}}
}
  \item{time}{
an optional vector of times in POSIXct format, where each time corresponds to a value in the \code{data} matrix or vector, used for pretty-printing in the \code{plot.oCPD} function.
}
  \item{hazard_func}{
hazard function used in the model. Defaults to a constant hazard, suitable for exponential family models.
}
  \item{m}{
initial value of \code{mu}, the mean of the data. Defaults to 0. As the mean is updated with every data point, this value does not need to be changed, but is safe to be experimented with.
}
  \item{k}{
initial value of \code{kappa}, basically a counter. Defaults to 0.01. May be useful to increase this to 1 if the data is large.
}
  \item{a}{
initial value of \code{alpha}, basically a counter. Defaults to 0.01. May be useful to increase this to 1 if the data is large.
}
  \item{b}{
initial value of \code{beta}, the variance of the data. Defaults to 1e-4. As the variance is updated with every data point, this value does not need to be changed, but is safe to be experimented with.
}
}
\details{
The primary result is a list of detected changepoints. Note that the list of changepoints must be interpreted; for example, sometimes the algorithm is unsure of the exact location of a change and prints several possible changepoints. Helper function \code{\link{findCP}}, called from \code{\link{plot.oCPD}}, will help reduce some of these changes.

The algorithm works by estimating the posterior distribution over the run-length, or the number of data points since the last changepoint. At each time, the run-length can either increase by one or reduce to zero.

The functions \code{\link{summary.oCPD}}, \code{\link{plot.oCPD}}, \code{\link{str.oCPD}} and \code{\link{print.oCPD}} are used to obtain summaries of the results.

See \code{\link{plot.oCPD}} for advice on how to interpret the results after plotting.
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
\seealso{\link{plot.oCPD}, \link{summary.oCPD}, \link{print.oCPD} for summaries of the results.}
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




