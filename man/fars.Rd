\name{fars}
\alias{fars}
\title{FARS Data}
\description{
  The \code{fars} data frame has 5662 rows and 5 columns of
  vehicle crashes in the United States in 2020.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{fips}{
      the five-digit codes of U.S. counties (FIPS 6-4).
    }
    \item{weekend}{
      a dummy variable for weekday (\code{0}) and weekend day (\code{1}).
    }
    \item{longitude}{
      longitude of the county's barycenter. 
    }
    \item{latitude}{
      latitude of the county's barycenter.
    }
	\item{count}{
	  the count of vehicle crashes.
	}
  }
}

\details{
The original data with individual accident records are publicly available from the Fatality Analysis Reporting System (FARS) database. This aggregated dataset is restricted to the year 2020, it excludes individual records with missing information on accident location. It also excludes accidents occurred in Alaska and Hawaii. The dataset is balanced, i.e., combinations of \code{fips} and \code{weekend} with no recorded accidents have 0 counts.
}

\source{
Fatality Analysis Reporting System, National Highway Traffic Safety Administration at nhtsa.gov.
}


\keyword{datasets}
