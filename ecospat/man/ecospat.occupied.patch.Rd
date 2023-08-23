\name{ecospat.occupied.patch}
\alias{ecospat.occupied.patch}

\title{
Extract occupied patches of a species in geographic space.)
}
\description{
This function determines the occupied patch of a species using standard IUCN criteria (AOO, EOO) or predictive binary maps from Species Distribution Models.
}
\usage{
ecospat.occupied.patch (bin.map, Sp.occ.xy, buffer = 0)
}
\arguments{ 
    \item{bin.map}{Binary map (single layer, raster stack or SpatRaster) from a Species Distribution Model.}
    \item{Sp.occ.xy}{xy-coordinates of the species presence.}
    \item{buffer}{numeric. Calculate occupied patch models from the binary map using a buffer (predicted area occupied by the species or within a buffer around the species, for details see ?extract).}
}

\details{
Predictive maps derived from SDMs inform about the potential distribution (or habitat suitability) of a species. Often it is useful to get information about the area of the potential distribution which is occupied by a species, e.g. for Red List assessments.
}
\value{
a RasterLayer or SpatRaster with value 1.
}
\author{
Frank Breiner \email{frank.breiner@wsl.ch} with the contribution of Flavien Collart
}
\references{

IUCN Standards and Petitions Subcommittee. 2016. Guidelines for Using the IUCN Red List Categories and Criteria. Version 12. Prepared by the Standards and Petitions Subcommittee. Downloadable from http://www.iucnredlist.org/documents/RedListGuidelines.pdf

}

\seealso{
   \code{\link{ecospat.rangesize}}, \code{\link{ecospat.mpa}}, \code{\link{ecospat.binary.model}}
}
\examples{

\donttest{
library(raster)
library(terra)
library(dismo)

# coordinates of the plots
xy <- ecospat.testData[,2:3]

# environmental data
predictors <- rast(get(load(system.file("extdata", "ecospat.testEnvRaster.Rdata", package="ecospat"))))
env <- terra::extract(predictors,xy,ID=F)  
spData <- cbind.data.frame(occ=ecospat.testData$Veronica_alpina,env)
mod <- glm(occ~ddeg0+I(ddeg0^2)+srad68+I(srad68^2),data=spData,family = binomial())

# predict to entire dataset
pred <- terra::predict(predictors,mod,type="response") 
plot(pred)
points(xy[spData$occ==1,])
  
### make binary maps
#arbitratry threshold
pred.bin.arbitrary <- ecospat.binary.model(pred,0.3)
names(pred.bin.arbitrary) <- "me.arbitrary"
# use MPA to convert suitability to binary map
mpa.cutoff <- ecospat.mpa2(pred,xy[spData$occ==1,])
pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)
names(pred.bin.mpa) <- "me.mpa"

mpa.ocp  <- ecospat.occupied.patch(pred.bin.mpa,xy[spData$occ==1,])
arbitrary.ocp  <- ecospat.occupied.patch(pred.bin.arbitrary,xy[spData$occ==1,])
par(mfrow=c(1,2))
plot(mpa.ocp) ## occupied patches: green area
points(xy[spData$occ==1,],col="red",cex=0.5,pch=19)
plot(arbitrary.ocp)
points(xy[spData$occ==1,],col="red",cex=0.5,pch=19)


## with buffer:
mpa.ocp  <- ecospat.occupied.patch(pred.bin.mpa,xy[spData$occ==1,], buffer=5000)
arbitrary.ocp  <- ecospat.occupied.patch(pred.bin.arbitrary,xy[spData$occ==1,], buffer=5000)

plot(mpa.ocp) ## occupied patches: green area
points(xy[spData$occ==1,],col="red",cex=0.5,pch=19)
plot(arbitrary.ocp)
points(xy[spData$occ==1,],col="red",cex=0.5,pch=19)
}
}
\keyword{file}
