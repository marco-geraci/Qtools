\name{predict.qrr}
\alias{predict.qrr}

\title{
Predictions from Quantile Ratio Regression Models
}

\description{
This function computes predictions based on quantile ratio regression models.
}

\usage{
\method{predict}{qrr}(object, newdata, na.action = na.pass,
type = "response", ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{qrr}.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
}

 \item{type}{
the type of prediction required. The default \code{"response"} is on the scale of the response variable, i.e. the values are back-transformed using the inverse of the link function \eqn{g^{-1}(Xb) = 1 + exp(Xb)}; the alternative \code{"link"} is on the scale of the linear predictor.
}
  \item{\dots}{
not used.
}
}


\value{
a vector of predictions.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{qrr}}
}

\keyword{ predict }

