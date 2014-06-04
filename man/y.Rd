\name{y}
\alias{y}
\docType{data}
\title{
Vector of observations
}
\description{
Vector of size 1600 which is a piecewise constant function corrupted by a Gaussian AR(1) process with rho=0.3 and where the standard deviation of the innovation process is equal to 0.1. The piecewise constant function has its breaks at the following positions: 222, 311, 711, 888, 1200, 1466 and takes the following values: 0, 1, 0, 1, 0, 1, 0 on each associated segment.
}
\usage{data(y)}
\format{
  The format is:
 num [1:1600] 0.11834 0.02428 0.00802 0.06716 0.10555 ...
}
\examples{
library(AR1seg)
data(y)
}
\keyword{datasets}
