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
xy, EOO, Model.within.eoo, AOO, resol, AOO.circles, d, lonlat, return.obj, save.obj, save.rangesize, directory)

ecospat.rangesize (bin.map = NULL,
                   ocp = T,
                   buffer = 0,
                   eoo.around.model = T,
                   eoo.around.modelocp = F,
                   xy = NULL,
                   EOO = T,
                   Model.within.eoo = T,
                   AOO = T,
                   resol = c(2000, 2000),
                   AOO.circles = F,
                   d = sqrt(2000^2/pi),
                   lonlat = FALSE,
                   return.obj = T,
                   save.obj = F,
                   save.rangesize = F,
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
