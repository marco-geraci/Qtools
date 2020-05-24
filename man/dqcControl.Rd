\name{dqcControl}
\alias{dqcControl}

\title{
Control parameters for dqc estimation
}
\description{
A list of parameters for controlling the fitting process.
}

\usage{
dqcControl(tau.range = c(0.001, 0.999), nt = 50, ndir = 50, seed = NULL)
}

\arguments{
	\item{tau.range}{
		vector with range of quantile probabilities.
	}

	\item{nt}{
		length of grid of quantiles within \code{tau.range}.
	}

	\item{ndir}{
		number of directions.
	}

	\item{seed}{
		seed for \code{\link{set.seed}}.
	}

}


\value{
a list of control parameters.
}

\references{
Viroli C, Farcomeni A, Geraci M (2020). Directional quantile-based classifiers (in preparation).
}

\author{
Marco Geraci}

\seealso{
\code{\link{dqc}}
}
\keyword{fitting}
\keyword{control}
