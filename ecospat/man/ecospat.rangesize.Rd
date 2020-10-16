\name{ecospat.rangesize}
\alias{ecospat.rangesize}

\title{
Quantification of the range size of a species using habitat suitability maps and IUCN criteria)
}
\description{
This function quantifies the range size of a species using standard IUCN criteria (Area of Occupancy AOO, Extent of Occurence EOO) or binary maps derived from Species Distribution Models.
}
\usage{
ecospat.rangesize (bin.map, ocp, buffer, eoo.around.model, eoo.around.modelocp, 
xy, EOO, Model.within.eoo, AOO, resol, AOO.circles, d, lonlat, return.obj, 
save.obj, save.rangesize, directory)

ecospat.rangesize (bin.map = NULL,
                   ocp = TRUE,
                   buffer = 0,
                   eoo.around.model = TRUE,
                   eoo.around.modelocp = FALSE,
                   xy = NULL,
                   EOO = TRUE,
                   Model.within.eoo = TRUE,
                   AOO = TRUE,
                   resol = c(2000, 2000),
                   AOO.circles = FALSE,
                   d = sqrt((2000 * 2)/pi),
                   lonlat = FALSE,
                   return.obj = TRUE,
                   save.obj = FALSE,
                   save.rangesize = FALSE,
                   directory = getwd())



}
\arguments{
    \item{bin.map}{Binary map (single layer or raster stack) from a Species Distribution Model.}
    \item{ocp}{logical. Calculate occupied patch models from the binary map (predicted area occupied by the species)}
    \item{buffer}{numeric. Calculate occupied patch models from the binary map using a buffer (predicted area occupied by the species or within a buffer around the species, for details see ?extract).}
    \item{eoo.around.model}{logical. The EOO around all positive predicted values from the binary map.}
    \item{eoo.around.modelocp}{logical. EOO around all positive predicted values of occupied patches.}
    \item{xy}{xy-coordinates of the species presence}
    \item{EOO}{logical. Extent of Occurrence. Convex Polygon around occurrences.}
    \item{Model.within.eoo}{logical. Area predicted as suitable by the model within EOO.}
    \item{AOO}{logical. Area of Occupancy ddervied by the occurrences.}
    \item{resol}{Resolution of the grid frame at which AOO should be calculated.}
    \item{AOO.circles}{logical. AOO calculated by circles around the occurrences instead of using a grid frame.}
    \item{d}{Radius of the AOO.circles around the occurrences.}
    \item{lonlat}{Are these longitude/latidue coordinates? (Default = FALSE).}
    \item{return.obj}{logical. should the objects created to estimate range size be returned (rasterfiles and spatial polygons). Default: TRUE}
    \item{save.obj}{logical. should objects be saved on hard drive?}
    \item{save.rangesize}{logical. should range size estimations be saved on hard drive .}
    \item{directory}{directory in which objects should be saved (Default = getwd())}
    }


\details{
The range size of a species is important for many conservation purposes, e.g. to assess the status of threat for IUCN Red Lists. This function quantifies the range size using different IUCN measures, i.e. the Area Of Occupancy (AOO), the Extent Of Occurrence (EOO) and from binary maps derived from Species Distribution Models (SDMs). Different ways to extract range size from SDMs are available, e.g. using occupied patches, the suitable habitat within EOO or a convex hull around the suitable habitat.
}
\value{
A list with the values of range size quantification and the stored objects used for quantification (of class RasterLayers, ahull, ConvexHull).
}
\author{
Frank Breiner \email{frank.breiner@wsl.ch}
}
\references{

IUCN. 2012. IUCN Red List Categories and Criteria: Version 3.1. Second edition. Gland, Switzerland and Cambridge, UK: IUCN. iv + 32pp.

IUCN Standards and Petitions Subcommittee. 2016. Guidelines for Using the IUCN Red List Categories and Criteria. Version 12. Prepared by the Standards and Petitions Subcommittee. Downloadable from http://www.iucnredlist.org/documents/RedListGuidelines.pdf

Pateiro-Lopez, B., and A. Rodriguez-Casal. 2010. Generalizing the Convex Hull of a Sample: The R Package alphahull. \emph{Journal of Statistical software}, \bold{34}, 1-28.
}

\seealso{
   \code{\link{ecospat.occupied.patch}}, \code{\link{ecospat.mpa}}, \code{\link{ecospat.binary.model}}
}

\keyword{file}
  
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

if (file.exists(path.to.renamed) & require(rJava)) {

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

### make binary maps

# use MPA to convert suitability to binary map
mpa.cutoff <- ecospat.mpa(pred,occ)

# use Boyce index to convert suitability to binary map
boyce <- ecospat.boyce(pred,  occ)
### use the boyce index to find a threshold
pred.bin.arbitrary <- ecospat.binary.model(pred,0.5)

pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)
names(pred.bin.mpa) <- "me.mpa"
pred.bin.arbitrary <- ecospat.binary.model(pred,0.5)
names(pred.bin.arbitrary) <- "me.arbitrary"

### rangesize calculations

library(rgeos)

rangesize <- ecospat.rangesize(stack(pred.bin.mpa,pred.bin.arbitrary),
                               xy=occ,
                               resol=c(1,1),
                               eoo.around.modelocp =TRUE,
                               AOO.circles = TRUE,
                               d=200000,
                               lonlat =TRUE)
rangesize$RangeSize

names(rangesize$RangeObjects)

par(mfrow=c(1,3))

plot(ecospat.binary.model(pred,0),legend=FALSE, main="IUCN criteria")

### IUCN criteria & derivates

# plot AOO
plot(rangesize$RangeObjects$AOO,add=TRUE, col="red",legend=FALSE)

# plot EOO
plot(rangesize$RangeObjects$EOO@polygons,add=TRUE, border="red", lwd=2)

# plot circles around occurrences
plot(rangesize$RangeObjects$AOO.circle@polygons,add=TRUE,border="blue")

for(i in 1:2){
  ## plot the occupied patches of the model
  plot(rangesize$RangeObjects$models.ocp[[i]],col=c("grey","blue","darkgreen"),
  main=names(rangesize$RangeObjects$models.ocp[[i]]),legend=FALSE)
  points(occ,col="red",cex=0.5,pch=19)
  ## plot EOO around model
  plot(rangesize$RangeObjects$eoo.around.model[[i]]@polygons,add=TRUE,border="blue",lwd=2)
  ## plot EOO around occupied patches
  plot(rangesize$RangeObjects$eoo.around.mo.ocp[[i]]@polygons,add=TRUE,border="darkgreen",
  lwd=2)
  ## plot the modeled area within EOO
  #plot(rangesize$RangeObjects$model.within.eoo[[i]],col=c("grey","blue","darkgreen"),legend=FALSE)
  #points(occ,col="red",cex=0.5,pch=19)
  #plot(rangesize$RangeObjects$EOO@polygons,add=TRUE, border="red", lwd=2)
}
par(mfrow=c(1,1))

### Alpha-hulls are not included in the function yet because of Licence limitations.
### However, alpha-hulls can easily be included manually (see also the help file of 
### the alpha hull package):

require(alphahull)
alpha = 2 # alpha value of 2 recommended by IUCN
  
del<-delvor(occ)
dv<-del$mesh
mn <- mean(sqrt(abs(del$mesh[,3]-del$mesh[,5])^2+abs(del$mesh[,4]-del$mesh[,6])^2))*alpha
alpha.hull<-ahull(del,alpha=mn) 
  
#Size of alpha-hulls
#areaahull(alpha.hull) #works but uses a deprecated function in alphahull 2.1

#plot alphahulls
plot(rangesize$RangeObjects$models.ocp[[i]],col=c("grey","blue","darkgreen"),
main=names(rangesize$RangeObjects$models.ocp[[i]]),legend=FALSE)
plot(alpha.hull,add=TRUE,lwd=1)
}
}
