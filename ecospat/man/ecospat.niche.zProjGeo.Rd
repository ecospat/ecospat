\name{ecospat.niche.zProjGeo}

\alias{ecospat.niche.zProjGeo}

\title{Projection of Occurrence Densities to the Geography}

\description{Creates a raster in geography with each pixel containing the occurrence densities extracted from a z object generated with \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.niche.zProjGeo(z1,env,cor)}

\arguments{
  \item{z1}{Species 1 occurrence density grid created by \code{ecospat.grid.clim.dyn}.}
  \item{env}{A SpatRaster of environmental variables corresponding to the background (\code{glob} in \code{ecospat.grid.clim.dyn}).}
  \item{cor}{FALSE by default. If TRUE corrects the occurrence densities of each species by the prevalence of the environments in their range}
}

\details{extracts the occurrence density of z objects created by \code{ecospat.grid.clim.dyn} for each point of the background (\code{glob}) using \code{extract} (package terra). The values are binded to the geographic coordinates of \code{env} and a raster is then recreated using \code{rast} with the argument type="xyz".
}

\value{raster of class RasterLayer}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}


\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}:481-497.


Petitpierre, B., C. Kueffer, O. Broennimann, C. Randin, C. Daehler and A. Guisan. 2012. Climatic niche shifts are rare among terrestrial plant invaders. \emph{Science}, \bold{335}:1344-1348.
}

\seealso{\code{\link{ecospat.plot.niche.dyn}}, \code{\link{ecospat.niche.dynIndexProjGeo}}}

\examples{
\donttest{

library(ade4)
library(terra)

data("ecospat.testNiche")
spp <- ecospat.testNiche
xy.sp1<-subset(spp,species=="sp1")[2:3] #Bromus_erectus

env<-rast(system.file("extdata","ecospat.testEnv.tif",package="ecospat"))
#plot(env)

env.sp1<-extract(env,xy.sp1)[,-1]
env.bkg<-na.exclude(values(env))

#################################### PCA-ENVIRONMENT ##################################

pca.cal <- dudi.pca(env.bkg, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# predict the scores on the axes
scores.bkg <- pca.cal$li  #scores for background climate
scores.sp1 <- suprow(pca.cal,env.sp1)$lisup				#scores for sp1

# calculation of occurence density (niche z)
  
z1 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg, scores.sp1,R=100)

plot(z1$z.uncor)
points(scores.sp1)

#################################### occurrence density in space ##################################

# sp1
geoz1<-ecospat.niche.zProjGeo(z1,env)
plot(geoz1,main="z1")
points(xy.sp1)
}}

