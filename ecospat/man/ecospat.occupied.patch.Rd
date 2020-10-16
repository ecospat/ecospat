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
    \item{bin.map}{Binary map (single layer or raster stack) from a Species Distribution Model.}
    \item{Sp.occ.xy}{xy-coordinates of the species presence.}
    \item{buffer}{numeric. Calculate occupied patch models from the binary map using a buffer (predicted area occupied by the species or within a buffer around the species, for details see ?extract).}
}

\details{
Predictive maps derived from SDMs inform about the potential distribution (or habitat suitability) of a species. Often it is useful to get information about the area of the potential distribution which is occupied by a species, e.g. for Red List assessments.
}
\value{
a RasterLayer with value 1.
}
\author{
Frank Breiner \email{frank.breiner@wsl.ch}
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
library(dismo)

### make a maxent model

# copy maxent.jar file in the right folder
path.from<-system.file("extdata", "maxent.txt", package="ecospat")
path.to <- paste0(system.file(package="dismo"), "/java/maxent.txt")
path.to.renamed <- paste0(system.file(package="dismo"), "/java/maxent.jar")
file.copy(path.from,path.to,overwrite = TRUE)
file.rename(path.to, path.to.renamed)

if (file.exists(path.to.renamed) & require(rJava) & require(igraph)) {

  # get predictor variables
  fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
                     pattern='grd', full.names=TRUE )
  predictors <- stack(fnames)
  #plot(predictors)

  # file with presence points
  occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
  occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
  colnames(occ) <- c("x","y")
  occ <- ecospat.occ.desaggregation(occ,min.dist=1)

  # fit a domain model, biome is a categorical variable
  me <- maxent(predictors, occ, factors='biome')

  # predict to entire dataset
  pred <- predict(me, predictors) 

  plot(pred)
  points(occ)
}

### to convert suitability to binary map

mpa.cutoff <- ecospat.mpa(pred,occ)

pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)
names(pred.bin.mpa) <- "me.mpa"
pred.bin.arbitrary <- ecospat.binary.model(pred,0.5)
names(pred.bin.arbitrary) <- "me.arbitrary"

### calculate occupied patch

mpa.ocp  <- ecospat.occupied.patch(pred.bin.mpa,occ)
arbitrary.ocp  <- ecospat.occupied.patch(pred.bin.arbitrary,occ)

par(mfrow=c(1,2))
plot(mpa.ocp) ## occupied patches: green area
points(occ,col="red",cex=0.5,pch=19)
plot(arbitrary.ocp)
points(occ,col="red",cex=0.5,pch=19)

## with buffer:
mpa.ocp  <- ecospat.occupied.patch(pred.bin.mpa,occ, buffer=500000)
arbitrary.ocp  <- ecospat.occupied.patch(pred.bin.arbitrary,occ, buffer=500000)

plot(mpa.ocp) ## occupied patches: green area
points(occ,col="red",cex=0.5,pch=19)
plot(arbitrary.ocp)
points(occ,col="red",cex=0.5,pch=19)
}

}
\keyword{file}
