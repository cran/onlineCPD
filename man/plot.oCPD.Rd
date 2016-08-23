\name{plot.oCPD}
\alias{plot.oCPD}
\title{
Plotting onlineCPD results
}
\description{
\code{plot.oCPD()} produces a summary plot of the results of \code{onlineCPD()}. Results are somewhat hard to interpret, see Details
}
\usage{
\method{plot}{oCPD}(x, lines = TRUE, title = "", leg.name = "Variable", cleanCP = TRUE, buffer = 10, ...)
}
\arguments{
  \item{x}{
the result of a call to \code{\link{onlineCPD}} or \code{\link{offlineCPD}}
}
  \item{lines}{
if true, plot red vertical lines in the top plot to denote detected changepoints
}
  \item{title}{
the title of the plot
}
  \item{leg.name}{
the title of the legend. Note that the names in the legend are taken from the column names in \code{x$data}.
}
  \item{cleanCP}{
if true, the function will call \code{\link{findCP}} to remove some changepoints
}
  \item{buffer}{
The minimum number of points that need to separate two changepoints before they are both replaced
}
  \item{\dots}{
(optional) additional arguments, ignored.
}
}
\details{
The plotted results can be difficult to interpret.

The top plot is the data plotted as a scatterplot, with each data column in a different colour. If \code{lines == TRUE} then red vertical lines will be plotted to delineate the different runs.

The bottom plot shows the algorithms results. The black smears show the probability that run length is a particular value on the y-axis at the time on the x-axis. That is, for each pair (x,y), the darkness is the probability that at time x the run length is y. If run continues, the smear continues to move diagonally up. If the run stops, the smear returns to zero on the y-axis. The red diagonal line plots the largest probability at each time. Note the log scale on Probability.

If \code{time} was not \code{NULL} in the call to \code{offlineCPD}, then \code{time} will be along the x-axis of both plots. 

See \code{\link{findCP}} for information on which changepoints will be removed if \code{cleanCP} is true.
}
\author{
Zachary Zanussi
}
\note{
Depends on ggplot2, reshape2, gridExtra and scales
}
\seealso{
\link{print.oCPD}, \link{summary.oCPD}, \link{str.oCPD} for summary results, and \code{\link{findCP}} for information on how changepoints are reduced
}
\examples{
set.seed(6)
x <- c(rnorm(50,mean=0.3,sd=0.15),rnorm(40,mean=0.7,sd=0.1),rnorm(60,mean=0.5,sd=0.15))
plot(offlineCPD(x))
}
