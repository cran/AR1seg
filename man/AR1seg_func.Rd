\name{AR1seg_func}
\alias{AR1seg_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmentation of an AR(1) Gaussian process
}
\description{
This function consists in an implementation of a robust approach to solve the problem
of multiple change-point estimation in the mean of a Gaussian AR(1) process. 
A robust estimator of the autoregression parameter is proposed and used to build a decorrelated series on which a classical penalized least-square approach is applied.
}
\usage{
AR1seg_func(y, Kmax = 15, rho = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Vector of observations
}
  \item{Kmax}{
Maximal number of segments
}
  \item{rho}{
It corresponds to the autoregression parameter. If it is equal to TRUE then it is estimated using a robust approach, otherwise the user has to provide a numerical value.
By default, the value of rho is TRUE.
}
}
\value{
Contains the following attributes:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{data }{Vector of observations}
\item{rho }{The estimator of rho if the argument rho=TRUE, otherwise the value provided by the user}
\item{decorrelated }{The decorrelated series using rho}
\item{breaks }{Matrix of size Kmax*Kmax. The line K=1,...,Kmax corresponds to the optimal segmentation of the series with K segments. By convention, the last break of each line is the length of the series.}
\item{selected }{Selected number of segments using the modified BIC criterion proposed by Zhang and Siegmund (2007)}
\item{SelectedBreaks }{Optimal segmentation with a number of segments equal to the value selected}
\item{PPbreaks }{Matrix of breaks obtained after the post-processing step}
\item{PPSelectedBreaks }{Result of the post-processing step applied to SelectedBreaks: it is the resulting segmentation of our approach}
\item{PPselected }{Length of the resulting segmentation (PPSelectedBreaks)}
\item{PPmean }{Empirical mean of the series on each segment of the resulting segmentation}
}
\references{
This function corresponds to the implementation of the robust approach for estimating 
change-points in the mean of an AR(1) Gaussian process by using the methodology described in the paper arXiv:1403.1958
}
\author{
S. Chakar, E. Lebarbier, C. Levy-Leduc, S. Robin
}
\note{
This package depends on the package Segmentor3IsBack
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, Kmax = 15, rho = TRUE) 
{
    l = length(y)
    if (rho) 
        rho = median((diff(y, lag = 2))^2)/median(diff(y)^2) - 
            1
    x = y[2:l] - rho * y[1:(l - 1)]
    S = Segmentor(x, model = 2, Kmax = Kmax)
    breaks = S@breaks
    for (i in 1:Kmax) {
        for (j in 1:i) breaks[i, j] = breaks[i, j] + 1
    }
    rm(i, j)
    parameters = S@parameters
    PP = function(t) {
        x = t
        l = length(x)
        i = 2
        while (l > 2 && i < l) {
            if (x[i] == x[i - 1] + 1 && x[i] != x[i + 1] - 1) {
                x = c(x[1:(i - 1)], x[(i + 1):l])
                l = l - 1
            }
            else i = i + 1
        }
        if (l > 1 && x[l - 1] == x[l] - 1) 
            x = x[1:(l - 1)]
        x
    }
    PPbreaks = matrix(0, nrow = Kmax, ncol = Kmax, dimnames = dimnames(breaks))
    PPbreaks[1, ] = breaks[1, ]
    for (i in 2:Kmax) {
        t = PP(breaks[i, 1:(i - 1)])
        PPbreaks[i, ] = c(t, l, rep(0, Kmax - length(t) - 1))
    }
    rm(i, t)
    fMa = function(t, mu) {
        M = c()
        t = c(0, t)
        for (i in 2:length(t)) {
            M = c(M, rep(mu[i - 1], t[i] - t[i - 1]))
        }
        M
    }
    sswg = function(br, param, series) {
        sum((series - fMa(br, param))^2)
    }
    sswgseg = function(seg, seri) {
        res = c()
        for (i in 1:(Kmax)) {
            res = c(res, sswg(seg@breaks[i, 1:i], seg@parameters[i, 
                1:i], seri))
        }
        res
    }
    minushalflogB = function(t, u) {
        t = t[t != 0]
        l = length(t)
        b = log(t[1]/u)
        if (l > 1) {
            for (i in 2:l) {
                b = b + log(t[i] - t[i - 1])
            }
        }
        b = -b/2
    }
    ZS = function(seg, seri) {
        u = length(seg@data)
        Kmax = seg@Kmax
        f = function(t) minushalflogB(t, u)
        wg = sswgseg(seg, seri)
        criterion = -(((u + 1):(u - Kmax + 2))/2) * log(wg) + 
            lgamma(((u + 1):(u - Kmax + 2))/2) - (0:(Kmax - 1)) * 
            log(u) + apply(seg@breaks, 1, f)
        selected = which.max(criterion)
        selected
    }
    selected = ZS(S, x)
    SelectedBreaks = breaks[selected, 1:selected]
    PPSelectedBreaks = PPbreaks[selected, ]
    PPSelectedBreaks = PPSelectedBreaks[PPSelectedBreaks != 0]
    PPselected = length(PPSelectedBreaks)
    vec1 = c(1, PPSelectedBreaks[1:(PPselected - 1)] + 1)
    vec2 = PPSelectedBreaks[1:(PPselected)]
    m = c()
    for (i in 1:PPselected) {
        m[i] = mean(y[vec1[i]:vec2[i]])
    }
    list(data = y, rho = rho, decorrelated = x, breaks = breaks, 
        PPbreaks = PPbreaks, selected = selected, SelectedBreaks = SelectedBreaks, 
        PPSelectedBreaks = PPSelectedBreaks, PPselected = PPselected, PPmean = m)
  }
}

