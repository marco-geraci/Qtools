\name{midq2q}
\alias{midq2q}

\title{
Recover Ordinary Conditional Quantiles from Conditional Mid-Quantiles}

\description{
This function recovers ordinary conditional quantile functions based on fitted mid-quantile regression models.
}

\usage{
midq2q(object, newdata, observed = FALSE, ...)
}

\arguments{
  \item{object}{
	an object of \code{\link{class}} \code{midrq}.
	}

  \item{newdata}{
	a required data frame in which to look for variables with which to predict.
	}
  \item{observed}{
	logical flag. If \code{TRUE}, ordinary quantiles are recovered from observed sample values. Otherwise, they are calcuated as rounded mid-quantiles. See details.
	}


 \item{\dots}{
not used.
}
}

\details{
If the values of the support of the random variable are equally spaced integers, then \code{observed} should ideally be set to \code{FALSE} so that the ordinary quantile is obtained by rounding the predicted mid-quantile. Otherwise, the function returns an integer observed in the sample. See Geraci and Farcomeni for more details.
}

\value{
a vector or a matrix of estimated ordinary quantiles. The attribute \code{Fhat} provides the corresponding estimated cumulative distribution.
}


\references{
Geraci, M. and A. Farcomeni. Mid-quantile regression for discrete responses. arXiv:1907.01945 [stat.ME]. URL: \url{https://arxiv.org/abs/1907.01945}.
}

\author{
Marco Geraci
}


\examples{

\dontrun{
# Esterase data
data(esterase)

# Fit quantiles 0.1, 0.15, ..., 0.85
fit <- midrq(Count ~ Esterase, tau = 2:17/20, data = esterase, type = 3, lambda = 0)

# Recover ordinary quantile function
xx <- seq(min(esterase$Esterase), max(esterase$Esterase), length = 5)
print(Qhat <- midq2q(fit, newdata = data.frame(Esterase = xx)))

# Plot
plot(Qhat, sub = TRUE)

}

}

\seealso{
\code{\link{plot.midq2q}}, \code{\link{predict.midrq}}
}

\keyword{predict}

