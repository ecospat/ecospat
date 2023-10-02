\name{ecospat.testEnv}

\docType{data}

\alias{ecospat.testEnv}

\title{Test Environmental SpatRasters for The Ecospat package}

\description{A stack of 5 topoclimatic SpatRasters at 250m resolution for the Western Swiss Alps. It includes "ddeg0" (growing degree-days above 0C), "mind68" (moisture index for month June to August),  "srad68" (solar radiation for month June to August), "slope25" (average of slopes at 25m resolution) and "topos25" (average of topographic positions at 25m resolution)}

\format{ecospat.testEnv is a tif file that contains the following SpatRasters:

 [1] "ddeg0"       
 [2] "mind68"      
 [3] "srad68"      
 [4] "slope25"            
 [5] "topos25"
}

\references{Zimmermann, N.E., F. Kienast. 2009. Predictive mapping of alpine grasslands in Switzerland: Species versus community approach. \emph{Journal of Vegetation Science}, \bold{10}, 469-482.
}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}

}

\examples{
\dontrun{
library(terra)
fpath <- system.file("extdata", "ecospat.testEnv.tif", package="ecospat")
env<-terra::rast(fpath)
plot(env)
}

}