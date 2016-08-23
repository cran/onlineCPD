\name{print.oCPD}
\alias{print.oCPD}
\alias{summary.oCPD}
\alias{str.oCPD}
\title{
Summarizing oCPD results
}
\description{
Summary, str and print methods for class \code{oCPD}
}
\usage{
\method{print}{oCPD}(x, ...)
\method{summary}{oCPD}(object, ...)
\method{str}{oCPD}(object, ...)
}
\arguments{
  \item{x}{
the result of a call to \code{\link{onlineCPD}} or \code{\link{offlineCPD}}
}
  \item{object}{ 
the result of a call to \code{\link{onlineCPD}} or \code{\link{offlineCPD}}
}
  \item{\dots}{
(optional) additional arguements, ignored.
}
}
\details{
The functions print and return invisibly the detected changepoints.
}
\author{
Zachary Zanussi
}



\seealso{
\code{\link{onlineCPD}} and \code{\link{plot.oCPD}}
}
\examples{
set.seed(6)
x <- c(rnorm(50,mean=0.3,sd=0.15),rnorm(40,mean=0.7,sd=0.1),rnorm(60,mean=0.5,sd=0.15))
res <- offlineCPD(x)
summary(res)
print(res)
str(res)
res
}
