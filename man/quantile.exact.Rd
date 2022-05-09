\name{quantile.exact}
\alias{quantile.exact}

\title{
Exact Confidence Intervals for Quantiles
}

\description{
Compute exact confidence intervals for quantiles using binomial probabilities
}

\usage{
quantile.exact(x, probs, level = 0.95)
}

\arguments{
  \item{x}{
numeric vector whose sample quantile and confidence intervals are to be calculated.
}
  \item{probs}{
numeric vector of probabilities with values in \code{[0,1]}.
}
  \item{level}{
nominal coverage level of the confidence interval.
}

}


\references{
Thompson W. R. On confidence ranges for the median and other expectation distributions for populations of unknown distribution form. The Annals of Mathematical Statistics 1936;7(3):122–128.
}

\author{
Marco Geraci
}

\examples{

x <- rnorm(100)
quantile.exact(x, p = c(0.1,0.5), level = 0.9)

}


\keyword{ unconditional quantiles }
