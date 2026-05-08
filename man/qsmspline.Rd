\name{qsmspline}
\alias{qsmspline}

\title{Quantile Regression with Smoothing Splines}
\description{
This function is used to fit a quantile  smoothing spline with L2 penalty.
}
\usage{
qsmspline(x, y, p, lambda = 1, maxIter = 300, eps = 1e-2, gamma = 10,
	aggressive = FALSE)
}

\arguments{
  \item{x,y}{
vectors giving the coordinates of the points to be interpolated.
}
  \item{p}{
quantile to be estimated.
}
  \item{lambda}{
penalty parameter.
}
  \item{maxIter}{
the maximum number of iterations.
}
  \item{eps}{
the absolute convergence tolerance. 
}
  \item{gamma}{
scaling for initial values.
}
  \item{aggressive}{
aggressive step size in Bosch et al (1995) - not yet implemented.  
}

}

\details{This is an implementation of Bosch et al's (1995) algorithm to fit a quantile smoothing spline with L2 penalty. The penalty parameter must be set by the user.
} 

\value{
a list containing the following components (see Bosch et al, 1995)

\item{a1}{first set of constraints.}
\item{a2}{second set of constraints.}
\item{b1}{Lagrangian associated with a1.}
\item{b2}{Lagrangian associated with a2.}
\item{fit}{fitted values.}
}


\references{
Bosch RJ, Ye Y, and Woodworth GG. A convergent algorithm for quantile regression with smoothing splines. Computational Statistics and Data Analysis. 1995;19(6):613-30.
}

\author{
Marco Geraci
}



\examples{
# Generate data
set.seed(123)
n <- 100
x <- sort(runif(n, 0, 2*pi))
y <- sin(x) + (1 + x)*rnorm(n, 0, 0.1)

# Fit median - lambda set at an arbitrary value
fit <- qsmspline(x = x, y = y, p = 0.5, lambda = 0.1)

# Plot
plot(x, y)
lines(x, fit[['fit']])
}

\keyword{conditional quantiles}
\keyword{smoothing}
