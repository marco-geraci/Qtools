\name{qrr}
\alias{qrr}

\title{
Quantile Ratio Regression
}
\description{
This function fits a quantile ratio regression model
}

\usage{
qrr(formula, data, taus, start = "rq", beta = NULL, 
tsf = "bc", symm = TRUE, dbounded = FALSE, linearize = TRUE, 
kernel = "Gaussian", maxIter = 10, epsilon = 1e-05, verbose = FALSE, 
method.rq = "fn", method.nlrq = "L-BFGS-B")
}

\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and the terms, separated by \code{+} operators, on the right.
}
  \item{data}{
a data frame in which to interpret the variables named in the formula.
}
  \item{taus}{
a vector of two quantiles for the ratio to be estimated (the order is irrelevant).
}
  \item{start}{
the algorithm with which obtain the starting values for one of the quantiles in the ratio. Possible options are \code{"rq"} (linear regression model -- see \code{\link[quantreg]{rq}}), \code{"tsrq"} (quantile regression transformation model -- see \code{\link{tsrq}}), \code{"conquer"} (fast linear regression model -- see \code{\link[conquer]{conquer}}), \code{"llqr"} (nonparametric linear regression model -- see \code{\link[quantdr]{llqr}})
}
  \item{beta}{
starting values for the regression coefficients. Set to 0 if left \code{NULL}.
}
  \item{tsf}{
if \code{start = "tsrq"}, see \code{\link{tsrq}}.
}
  \item{symm}{
if \code{start = "tsrq"}, see \code{\link{tsrq}}.
}
  \item{dbounded}{
if \code{start = "tsrq"}, see \code{\link{tsrq}}.
}
  \item{linearize}{
logical flag. If \code{TRUE} (default), estimation is carried out with the linearized iterative algorithm of Farcomeni and Geraci (2023) by repeated calls to an appropriate linear estimation algorithm. Otherwise, the algorithm calls a nonlinear estimation routine. See argument \code{method.rq} and \code{method.nlrq} further below.
}
  \item{kernel}{
an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. 
}
  \item{maxIter}{
maximum number of iterations for fitting.
}
  \item{epsilon}{
tolerance for convergence.
}
  \item{verbose}{
logical flag. If \code{TRUE}, progress on estimation is print out.
}
  \item{method.rq}{
the method used to compute the linear fit. If \code{linearize = TRUE}, the options are \code{"conquer"} or any of those from \code{\link[quantreg]{rq}} (see the argument \code{method}).
}
  \item{method.nlrq}{
the method used to compute the nonlinear fit. If \code{linearize = FALSE}, the options are those from \code{\link[quantreg]{nlrq}} (see the argument \code{method}).
}

}

\details{
These function implements quantile ratio regression as discussed by Farcomeni and Geraci (see references). The general model is assumed to be \eqn{g(Q_{Y|X}(\tau_{1})/Q_{Y|X}(\tau_{2})) = \eta = Xb} where \eqn{Q} denotes the conditional quantile function, \eqn{Y} is the response variable, \eqn{X} a design matrix, \eqn{g} is a monotone link function, and \eqn{\tau_{1}} and \eqn{\tau_{2}} the levels of the two quantiles in the ratio. In the current implementation, \eqn{g(u) = log(u - 1)}, which ensures monotonocity (non-crossing) of the quantiles and leads to the familiar interpretation of the inverse logistic transformation. 
}

\references{
Farcomeni A. and Geraci M. Quantile ratio regression. 2023. Working Paper.
}

\author{
Marco Geraci
}

\examples{

set.seed(123)
n <- 5000
x <- runif(n, -0.5, 0.5)
R <- 1 + exp(0.5 + 0.5*x)

# fit quintile ratio regression
alpha <- 1/log(R)*log(log(1-0.8)/log(1-0.2))
y <- rweibull(n, shape = alpha, scale = 1)
dd <- data.frame(x = x, y = y)
qrr(y ~ x, data = dd, taus = c(.2,.8))

# fit Palma ratio regression
alpha <- 1/log(R)*log(log(1-0.9)/log(1-0.4))
y <- rweibull(n, shape = alpha, scale = 1)
dd <- data.frame(x = x, y = y)
qrr(y ~ x, data = dd, taus = c(.4,.9))

}


\seealso{
\code{\link{coef.qrr}}, \code{\link{predict.qrr}}, \code{\link{summary.qrr}}, \code{\link{vcov.qrr}}
}

\keyword{ quantile ratios }
\keyword{ conditional quantiles }

