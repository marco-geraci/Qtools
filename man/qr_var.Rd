\name{qr_var}
\alias{qtr_var}
\alias{qr_var}

\title{
Variance for Ratios of Quantiles and Tail means
}
\description{
This functions returns the variance of the ratio between quantiles and between tail means.
}

\usage{
qr_var(object, newdata, se = c("iid", "nid"))
qtr_var(fitL, fitU, newdata, se = c("iid", "nid"),
wL = NULL, wU = NULL)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{rq} with exactly two fitted regression quantiles.
}

  \item{newdata}{
an optional data frame in which to look for variables with which to predict.
}

  \item{se}{
the standard error method for individual regression quantiles. Only options allowed are \code{iid} and \code{nid}. See \code{\link{summary.rq}}.
}

  \item{fitL}{
an object of \code{\link{class}} \code{rq} with a grid of fitted regression quantiles for the lower tail.
}

  \item{fitU}{
an object of \code{\link{class}} \code{rq} with a grid of fitted regression quantiles for the upper tail.
}

  \item{wL}{
weights for the lower tail used for averaging. If not provided, uniform weights are used.
}

  \item{wU}{
weights for the upper tail used for averaging. If not provided, uniform weights are used.
}

}

\details{
The function \code{qr_var} gives the large-sample approximation of the standard error of the ratio of two predicted quantiles. Similarly, the function \code{qtr_var} gives the variance of the L-statistic defined by the ratio of tail integrals. See examples for further details. 
}

\references{
Geraci M. Applied Statistics Using Quantiles. To appear.
}


\author{
Marco Geraci
}

\examples{

library(quantreg)
## Ratio of two quantiles
set.seed(123)
n <- 5000
x <- runif(n)
e <- rnorm(n)
y <- 5 + 2*x + (1 + x)*e
newd <- data.frame(x = 0.5)

fit <- rq(y ~ x, tau = c(0.1, 0.9), method = "fn")
Q <- predict(fit, newd)

# ratio of quantiles
QR <- Q[2]/Q[1]
QR
# standard error
v <- qr_var(fit, newdata = newd, se = "nid")
sqrt(v)

## Ratio of two tail means

set.seed(123)
n <- 5000
pL <- seq(0.01, 0.10, by = 0.01)
pU <- seq(0.90, 0.99, by = 0.01)
dpL <- diff(c(0, pL))
dpU <- diff(c(pU, 1))

x <- runif(n)
e <- rnorm(n)
y <- 5 + 2*x + (1 + x)*e
newd <- data.frame(x = 0.5)

fitL <- rq(y ~ x, tau = pL, method = "fn")
fitU <- rq(y ~ x, tau = pU, method = "fn")

QL <- sum(predict(fitL, newd)*dpL)/sum(dpL)
QU <- sum(predict(fitU, newd)*dpU)/sum(dpU)

# quantile tail ratio
QTR <- QU/QL
QTR
# standard error
v <- qtr_var(fitL, fitU, newdata = newd, se = "nid", wL = dpL, wU = dpU)
sqrt(v)

}

\seealso{
\code{\link{qrr}}
}

\keyword{ summary }
