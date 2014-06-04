\name{AR1seg-package}
\alias{AR1seg-package}
\alias{AR1seg}
\docType{package}
\title{
Segmentation of an AR(1) Gaussian process
}
\description{
This package consists in an implementation of a robust approach to solve the problem
of multiple change-point estimation in the mean of a Gaussian AR(1) process. 
A robust estimator of the autoregression parameter is proposed and used to build a decorrelated series on which a classical penalized least-square approach is applied.
}
\details{
\tabular{ll}{
Package: \tab AR1seg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-04\cr
License: \tab GPL-2\cr
}
}
\author{
S. Chakar, E. Lebarbier, C. Levy-Leduc, S. Robin

Maintainer: Souhil Chakar <souhil.chakar@agroparistech.fr>
}
\references{
S. Chakar, E. Lebarbier, C. Levy-Leduc, S. Robin. A robust approach to multiple change-point estimation in an AR(1) process, arXiv:1403.1958. 
}

\examples{
library(AR1seg)
data(y)
res=AR1seg_func(y,Kmax=15,rho=TRUE)
a=c(1,res$PPSelectedBreaks[1:(res$PPselected-1)]+1)
b=res$PPSelectedBreaks[1:(res$PPselected)]
Bounds=cbind(a,b)
mu.fit=rep(res$PPmean,Bounds[,2]-Bounds[,1]+1)
plot(y)
lines(mu.fit,col="red")
}
