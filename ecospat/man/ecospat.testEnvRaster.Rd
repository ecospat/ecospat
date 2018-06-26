\name{ecospat.testEnvRaster}

\docType{data}

\alias{ecospat.testEnvRaster}

\title{Test Environmental Rasters for The Ecospat package}

\description{A stack of 5 topoclimatic rasters at 250m resolution for the Western Swiss Alps. It includes "ddeg0" (growing degree-days above 0C), "mind68" (moisture index for month June to August),  "srad68" (solar radiation for month June to August), "slope25" (average of slopes at 25m resolution) and "topos25" (average of topographic positions at 25m resolution)}

\format{ecospat.testEnvRaster is a RasterBrick encapsulated in a .Rdata that contains the following rasters:

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
fpath <- system.file("extdata", "ecospat.testEnvRaster.RData", package="ecospat")
load(fpath)
plot(env)
}

}