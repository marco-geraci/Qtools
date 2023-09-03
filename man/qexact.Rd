\name{qexact}
\alias{qexact}

\title{
Exact Confidence Intervals for Quantiles
}

\description{
Compute exact confidence intervals for quantiles of continuous random variables using binomial probabilities
}

\usage{
qexact(x, probs = 0.5, level = 0.95)
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

\details{

This function calculates exact confidence intervals for quantiles at level \code{probs} from a vector \code{x} of length \code{n}. It does so by first determining the confidence level for all possible pairwise combinations of order statistics from 1 to \code{n}. This entails "\code{n} choose \code{2}" possible confidence intervals before selecting the one with the level closest to \code{level}. If the procedure yields more than one such confidence intervals, then the interval with smallest width is returned.

Caution: for large \code{n}, the procedure may reach the limit on the number of nested expressions. See \code{gtools::combinations} and \code{\link{options}(expressions)} for additional information. However, if you have a large \code{n}, then consider estimating an asymptotic approximation of the confidence interval.
}


\references{
Thompson W. R. On confidence ranges for the median and other expectation distributions for populations of unknown distribution form. The Annals of Mathematical Statistics 1936;7(3):122-128.
}


\author{
Marco Geraci
}

\examples{
x <- rnorm(100)
qexact(x, p = c(0.1,0.5), level = 0.9)
}


\keyword{ unconditional quantiles }
