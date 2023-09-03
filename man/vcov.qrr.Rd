\name{vcov.qrr}
\alias{vcov.qrr}

\title{
Variance-Covariance Matrix for a Fitted Quantile Ratio Regression Model Object
}
\description{
This functions returns the variance-covariance matrix of the coefficients of a fitted \code{qrr} model object.
}

\usage{
\method{vcov}{qrr}(object, method = "approximate", R = 200, update = TRUE, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{qrr}.
}
  \item{method}{
if \code{"approximate"}, the variance-covariance estimate is approximated by the inverse of the numerical Hessian. The latter is calculated as detailed in Farcomeni and Geraci (2023). If \code{"boot"}, the variance-covariance estimate is calculated by means of ordinary bootstrap (see \code{\link[boot]{boot}}).
}
  \item{R}{
the number of bootstrap replications.
}
  \item{update}{
logical flag. If \code{TRUE} (the default), the statistic to be resampled is obtained via an \code{\link{update}} of the \code{qrr} object. If \code{FALSE}, then the statistic to be resampled is obtained via a \code{\link{do.call}} of the \code{qrr} object. See details.
}

  \item{...}{
not used.
}

}

\details{
The use of \code{update = FALSE} is preferred when the function \code{vcov.qrr} is called from within another function.
}

\references{
Farcomeni A. and Geraci M. Quantile ratio regression. 2023. Working Paper.
}


\author{
Marco Geraci with contributions from Alessio Farcomeni
}

\seealso{
\code{\link{qrr}}
}

\keyword{ summary }
