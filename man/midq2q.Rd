\name{midq2q.midrq}
\alias{midq2q.midrq}
\alias{midq2q.midquantile}
\alias{midq2q}

\title{
Recover Ordinary Quantiles from Mid-Quantiles}

\description{
This function recovers ordinary quantiles from fitted mid-quantile objects.
}

\usage{
\method{midq2q}{midquantile}(object, observed = FALSE, ...)
\method{midq2q}{midrq}(object, observed = FALSE, ..., newdata, offset, na.action = na.pass)
}

\arguments{
  \item{object}{
	an object of \code{\link{class}} \code{midquantile} or \code{midrq}.
	}
  \item{observed}{
	logical flag. If \code{TRUE}, ordinary quantiles are recovered from observed sample values. Otherwise, they are calcuated as rounded mid-quantiles. See details.
	}

  \item{newdata}{
	optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
	}
  \item{offset}{
	an optional offset to be included in the model frame (when \code{newdata} is provided).
	}
  \item{na.action}{
	function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
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

# Fit conditional mid-quantiles 0.1, 0.15, ..., 0.85
fit <- midquantile(esterase$Count, probs = 2:17/20)

# Recover ordinary quantile function
print(Qhat <- midq2q(fit))

# Plot
plot(Qhat)

# Fit conditional mid-quantiles 0.1, 0.15, ..., 0.85
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

