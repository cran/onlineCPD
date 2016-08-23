\name{studentpdf}
\alias{studentpdf}
\title{
Student Probability Density Function
}
\description{
Calculate Student PDF; for use inside \code{\link{onlineCPD}} and \code{\link{offlineCPD}}

}
\usage{
studentpdf(x, mu, var, nu)
}
\arguments{
  \item{x}{
data point
}
  \item{mu}{
value of mu, the mean
}
  \item{var}{
the variance
}
  \item{nu}{
value of nu, the number of degrees of freedom
}
}
\details{
Used in the \code{\link{onlineCPD}} and \code{\link{offlineCPD}} function
}
\value{
The sample of Student's t distribution
}
\author{
Zachary Zanussi
}
