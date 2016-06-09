\name{ecospat.sample.envar}

\alias{ecospat.sample.envar}

\title{Sample Environmental Variables}

\description{Add environmental values to a species dataframe.}

\usage{ecospat.sample.envar (dfsp, colspxy, colspkept = "xy", dfvar, 
colvarxy, colvar = "all", resolution)}

\arguments{
  \item{dfsp}{A species dataframe with x (long), y (lat) and optional other variables.}
  \item{colspxy}{The range of columns for x (long) and y (lat) in dfsp.}
  \item{colspkept}{The columns of dfsp that should be kept in the final dataframe (default: xy).}
  \item{dfvar}{A dataframe object with x, y and environmental variables.}
  \item{colvarxy}{The range of columns for x and y in dfvar.}
  \item{colvar}{The range of enviromental variable columns in dfvar (default: all except xy).}
  \item{resolution}{The distance between x,y of species and environmental datafreme beyond which values shouldn't be added.} 
}

\details{The xy (lat/long) coordinates of the species occurrences are compared to those of the environment dataframe and the value of the closest pixel is added to the species dataframe. When the closest environment pixel is more distant than the given resolution, NA is added instead of the value. This function is similar to sample() in ArcGIS.}

\value{A Dataframe with the same rows as dfsp, with environmental values from dfvar in column.}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}


\examples{
\dontrun{
spp <- ecospat.testNiche
sp1 <- spp[1:32,1:3]
occ.sp1 <- ecospat.occ.desaggregation(dfvar=sp1,colxy=2:3,colvar=NULL, min.dist=500,plot=TRUE)
clim <- ecospat.testData[2:8]

occ_sp1 <- na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2,colspkept=1:2,dfvar=clim,
colvarxy=1:2,colvar="all",resolution=25))
}
}