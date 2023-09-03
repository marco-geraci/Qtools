\name{summary.qrr}
\alias{summary.qrr}

\title{
Summary for Quantile Ratio Regression Models
}
\description{
This functions gives a summary list for a quantile ratio regression model.
}

\usage{
\method{summary}{qrr}(object, se = "approximate", R = 200,
update = TRUE, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{summary.qrr}.
}
  \item{se}{
specifies the method used to compute standard errors. See argument \code{method} in \code{\link{vcov.qrr}}.
}
  \item{R}{
number of bootstrap replications.
}
  \item{update}{
see argument \code{update} in \code{\link{vcov.qrr}}.
}
  \item{\dots}{
not used.
}
}

\references{
Farcomeni A. and Geraci M. Quantile ratio regression. 2023. Working Paper.
}


\author{
Marco Geraci
}

\seealso{
\code{\link{qrr}}
}


\keyword{ bootstrap }
\keyword{ summary }
