\name{ecospat.testNiche.inv}

\docType{data}

\alias{ecospat.testNiche.inv}

\title{Test Data For The Niche Dynamics Analysis In The Invaded Range Of A Hypothetical Species}

\description{Data frame that contains geographical coordinates, environmental variables, occurrence sites for the studied species and the prediction of its distribution in the invaded range. These predictions are provided by SDM calibrated on the native range.}

\usage{data(ecospat.testNiche.inv)}

\format{ecospat.testNiche.inv is a data frame with the following columns:
  \describe{
    \item{\code{x}}{Longitude, in WGS84 coordinate system of the species occurrence.}
    \item{\code{y}}{Latitude, in WGS84 coordinate system of the species occurrence.}
    \item{\code{aetpet}}{Ratio of actual to potential evapotranspiration.}
    \item{\code{gdd}}{Growing degree-days above 5 degrees C.}
    \item{\code{p}}{Annual amount of precipitations.}
    \item{\code{pet}}{Potential evapotranspiration.}
    \item{\code{stdp}}{Annual variation of precipitations.}
    \item{\code{tmax}}{Maximum temperature of the warmest month.}
    \item{\code{tmin}}{Minimum temperature of the coldest month.}
    \item{\code{tmp}}{Annual mean temperature.}
    \item{\code{species_occ}}{Presence records of the species occurrence.}
    \item{\code{predictions}}{Species Distribution Model predictions of the studied species.}
}
}
\details{
The study area is Australia, which is the invaded range of the hypothetical species.

Eight topo-climatic explanatory variables to quantify niche differences: ratio of the actual potential evapotranspiration; growing degree days; precipitation; potential evapotranspiration; annual variation of precipitations; maximum temperature of the warmest month; minimum temperature of the coldest month; and annual mean temperature.
}

\references{Petitpierre, B., C. Kueffer, O. Broennimann, C. Randin, C. Daehler and A. Guisan. 2012. Climatic niche shifts are rare among terrestrial plant invaders. \emph{Science}, \bold{335}, 1344-1348.}


\author{Blaise Petitpierre \email{bpetitpierre@gmail.com} and Valeria Di Cola \email{valeria.dicola@unil.ch}
}

\seealso{
   \code{\link{ecospat.testNiche.nat}}
}

\examples{
data(ecospat.testNiche.inv)
str(ecospat.testNiche.inv)
dim(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)
}
