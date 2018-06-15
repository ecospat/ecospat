\name{ecospat.testMdr}

\docType{data}

\alias{ecospat.testMdr}

\title{Test Data For The ecospat.mdr function}

\description{Data frame that contains presence records the species \code{Centaurea stoebe} along years in North America.
}

\usage{data("ecospat.testMdr")}

\format{A data frame with 102 observations of \code{Centaurea stoebe}.
  \describe{
    \item{\code{latitude}}{Latitude, in WGS coordinate system.}
    \item{\code{longitude}}{Longitude, in WGS coordinate system.}
    \item{\code{date}}{Year of the presence record.}
}
}
\details{
Simplified dataset to exemplify the use of the ecospat.mdr function to calculate minimum dispersal routes.
}


\references{
Broennimann, O., P. Mraz, B. Petitpierre, A. Guisan, and H. Muller-Scharer. 2014. Contrasting spatio-temporal climatic niche dynamics during the eastern and western invasions of spotted knapweed in North America.\emph{Journal of biogeography}, \bold{41}, 1126-1136.

Hordijk, W. and O. Broennimann. 2012. Dispersal routes reconstruction and the minimum cost arborescence problem. \emph{Journal of theoretical biology}, \bold{308}, 115-122.
}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}
}

\examples{
data(ecospat.testMdr)
str(ecospat.testMdr)
dim(ecospat.testMdr)
}
