\name{ecospat.niche.dynIndexProjGeo}

\alias{ecospat.niche.dynIndexProjGeo}

\title{Projection of niche dynamic indices to the Geography}

\description{Creates a SpatRaster in geography with each pixel containing a
  niche dynamic index (stability, expansion, or unfilling) extracted
  from 2 niches generated with \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.niche.dynIndexProjGeo(z1, z2, proj, env)}

\arguments{
  \item{z1}{Species 1 occurrence density grid created
    by \code{ecospat.grid.clim.dyn}.}
  \item{z2}{Species 2 occurrence density grid created
    by \code{ecospat.grid.clim.dyn}.}
  \item{proj}{type of projection (0,1 or 2)}
  \item{env}{A SpatRaster of environmental variables
    corresponding to the background (\code{glob}
    in \code{ecospat.grid.clim.dyn}).}
}

\details{Extracts the niche dynamic index of objects created
  by \code{ecospat.niche.dyn.index} for each background point
  (\code{glob}) using \code{extract} from the terra package. The values
  are projected to the geographic coordinates of \code{env} and returned
  as a SpatRaster layer. If proj=0, the common background for z1 and z2 (\code{glob})
  is used. If proj=1, the background (\code{glob1}) of z1 is used. If proj=2, 
  the background (\code{glob1}) of z2 is used.
}

\value{a SpatRaster with cell values 0-3. 1: unfilling
  (habitat occupied by species 1); 2: expansion (habitat occupied only
  by species 2); 3: stability (habitat occupied by both species); and 0
  (unoccupied habitat).}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}, Tyler
  Smith \email{tyler@plantarum.ca}}

\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B.
  Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C.
  Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring
  ecological niche overlap from occurrence and spatial environmental
  data. \emph{Global Ecology and Biogeography}, \bold{21}:481-497.


  Petitpierre, B., C. Kueffer, O. Broennimann, C. Randin, C. Daehler and
  A. Guisan. 2012. Climatic niche shifts are rare among terrestrial
  plant invaders. \emph{Science}, \bold{335}:1344-1348.
}

\seealso{\code{\link{ecospat.plot.niche.dyn}},
  \code{\link{ecospat.niche.dyn.index}},
  \code{\link{ecospat.niche.zProjGeo}}}

\examples{
\donttest{

library(ade4)
library(terra)

data(ecospat.testNiche)
spp <- ecospat.testNiche
xy.sp1 <- subset(spp, species=="sp1")[2:3] #Bromus_erectus
xy.sp2 <- subset(spp, species=="sp4")[2:3] #Pritzelago_alpina

env<-terra::rast(system.file("extdata","ecospat.testEnv.tif",package="ecospat"))

env.sp1 <- terra::extract(env, xy.sp1)[,-1]
env.sp2 <- terra::extract(env, xy.sp2)[,-1]
env.bkg <- na.exclude(terra::values(env))

#####################
## PCA-ENVIRONMENT ##
#####################

library(ade4)
pca.cal <- ade4::dudi.pca(env.bkg, center = TRUE, scale = TRUE,
                    scannf = FALSE, nf = 2)

# predict the scores on the axes
scores.bkg <- pca.cal$li  #scores for background climate
scores.sp1 <- ade4::suprow(pca.cal,env.sp1)$lisup  #scores for sp1
scores.sp2 <- ade4::suprow(pca.cal,env.sp2)$lisup  #scores for sp2

# calculation of occurence density (niche z)
  
z1 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg,
                            scores.sp1, R = 100)
z2 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg,
                            scores.sp2, R = 100)

plot(z1$z.uncor)
points(scores.sp1, pch = 21, bg = "green")

plot(z2$z.uncor)
points(scores.sp2, pch = 21, bg = "red")

ecospat.plot.niche.dyn(z1, z2)
points(scores.sp1, pch = 21, bg = "green")
points(scores.sp2, pch = 21, bg = "red")
ecospat.niche.overlap(z1,z2 ,cor = TRUE)

##############################
## dynamic indices in space ##
##############################

geozS <- ecospat.niche.dynIndexProjGeo(z1, z2, proj=0, env)

plot(geozS, col = c("grey", "green", "red", "blue"),
     legend = FALSE)
points(xy.sp1, bg = "green", pch = 21)
points(xy.sp2, bg = "red", pch = 21)
}}
