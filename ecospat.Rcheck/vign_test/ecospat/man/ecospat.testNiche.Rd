\name{ecospat.testNiche}

\docType{data}

\alias{ecospat.testNiche}

\title{Test Data For The Niche Overlap Analysis}

\description{Data frame that contains occurrence sites for each species, long, lat and the name of the species at each site.}

\usage{data(ecospat.testNiche)}

\format{ecospat.testNiche is a data frame with the following columns:
  \describe{
    \item{\code{species}}{sp1, sp2, sp3 and sp4.}
    \item{\code{long}}{Longitude, in Swiss plane coordinate system of the vegetation plot.}
    \item{\code{lat}}{Latitude, in Swiss plane coordinate system of the vegetation plot.}
    \item{\code{Spp}}{Scientific name of the species used in the exmaple: Bromus_erectus_sstr, Saxifraga_oppositifolia, Daucus_carota and Pritzelago_alpina_sstr.}
}
}
\details{
List of occurence sites for the species.
}

\author{Antoine Guisan \email{antoine.guisan@unil.ch}, Anne Dubuis \email{anne.dubuis@gmail.com} and Valeria Di Cola \email{valeria.dicola@unil.ch}
}

\seealso{
   \code{\link{ecospat.testData}}
}

\examples{
data(ecospat.testNiche)
dim(ecospat.testNiche)
names(ecospat.testNiche)
}
