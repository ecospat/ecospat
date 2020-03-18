\name{ecospat.grid.clim.dyn}

\alias{ecospat.grid.clim.dyn}

\title{Dynamic Occurrence Densities Grid}

\description{Create a grid with occurrence densities along one or two gridded environmental gradients.}

\usage{ecospat.grid.clim.dyn (glob, glob1, sp, R, th.sp, th.env, geomask)}

\arguments{
  \item{glob}{A two-column dataframe (or a vector) of the environmental values (in column) for background pixels of the whole study area (in row).}
  \item{glob1}{A two-column dataframe (or a vector) of the environmental values (in column) for the background pixels of the species (in row).}
  \item{sp}{A two-column dataframe (or a vector) of the environmental values (in column) for the occurrences of the species (in row).}
  \item{R}{The resolution of the grid.}
  \item{th.sp}{The quantile used to delimit a threshold to exclude low species density values.}
  \item{th.env}{The quantile used to delimit a threshold to exclude low environmental density values of the study area.}
  \item{geomask}{A geographical mask to delimit the background extent if the analysis takes place in the geographical space.}

 It can be a SpatialPolygon or a raster object. Note that the CRS should be the same as the one used for the points.


}

\details{Using the scores of an ordination (or SDM prediction), create a grid z of RxR pixels (or a vector of R pixels when using scores of dimension 1 or SDM predictions) with occurrence densities. Only scores of one, or two dimensions can be used.
\code{th.sp} is the quantile of the distribution of species density at occurrence sites. 
For example, if \code{th.sp} is set to 0.05, the the species niche is drawn by including 95 percent of the species occurrences, removing the more marginal populations. 
Similarly, \code{th.env} is the quantile of the distribution of the environmental density at all sites of the study area. 
If \code{th.env} is set to 0.05, the delineation of the study area in the environmental space includes 95 percent  of the study area, removing the more marginal sites of the study area. 
By default, these thresholds are set to 0 but can be modified, depending on the importance of some marginal sites in the delineation of the species niche and/or the study area in the environmnental space. It is recommended to check if the shape of the delineated niche and study area corresponds to the shape of the plot of the PCA scores (or any other ordination techniques used to set the environmental space). 
Visualisation of the gridded environmental space can be done through the functions \code{\link[ecospat]{ecospat.plot.niche}} or \code{\link[ecospat]{ecospat.plot.niche.dyn}}
If you encounter a problem during your analyses, please first read the FAQ section of "Niche overlap" in http://www.unil.ch/ecospat/home/menuguid/ecospat-resources/tools.html
The argument \code{geomask} can be a SpatialPolygon or a raster object. 
Note that the CRS should be the same as the one used for the points.
}

\value{A grid z of RxR pixels (or a vector of R pixels) with z.uncor being the density of occurrence of the species, and z.cor the occupancy of the environment by the species (density of occurrences divided by the desinty of environment in the study area.}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch} and Blaise Petitpierre \email{bpetitpierre@gmail.com}}


\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}:481-497.


Petitpierre, B., C. Kueffer, O. Broennimann, C. Randin, C. Daehler and A. Guisan. 2012. Climatic niche shifts are rare among terrestrial plant invaders. \emph{Science}, \bold{335}:1344-1348.
}


\seealso{\code{\link{ecospat.plot.niche.dyn}}}

\examples{
\donttest{

data(ecospat.testNiche)
data(ecospat.testData)

spp <- ecospat.testNiche
clim <- ecospat.testData[2:8]

occ.sp_test <- na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=2:3,colspkept=1:3,dfvar=clim,
colvarxy=1:2,colvar="all",resolution=25))

occ.sp<-cbind(occ.sp_test,spp[,4]) #add species names

# list of species
sp.list<-levels(occ.sp[,1])
sp.nbocc<-c()

for (i in 1:length(sp.list)){sp.nbocc<-c(sp.nbocc,length(which(occ.sp[,1] == sp.list[i])))} 
#calculate the nb of occurences per species

sp.list <- sp.list[sp.nbocc>4] # remove species with less than 5 occurences
nb.sp <- length(sp.list) #nb of species
ls()
# selection of variables to include in the analyses 
# try with all and then try only worldclim Variables
Xvar <- c(3:7)
nvar <- length(Xvar)

#number of interation for the tests of equivalency and similarity
iterations <- 100
#resolution of the gridding of the climate space
R <- 100
#################################### PCA-ENVIRONMENT ##################################
data<-rbind(occ.sp[,Xvar+1],clim[,Xvar]) 
w <- c(rep(0,nrow(occ.sp)),rep(1,nrow(clim)))
pca.cal <- dudi.pca(data, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

####### selection of species ######
sp.list
sp.combn <- combn(1:2,2)

for(i in 1:ncol(sp.combn)) {
  row.sp1 <- which(occ.sp[,1] == sp.list[sp.combn[1,i]]) # rows in data corresponding to sp1
  row.sp2 <- which(occ.sp[,1] == sp.list[sp.combn[2,i]]) # rows in data corresponding to sp2
  name.sp1 <- sp.list[sp.combn[1,i]]
  name.sp2 <- sp.list[sp.combn[2,i]]
  # predict the scores on the axes
  scores.clim <- pca.cal$li[(nrow(occ.sp)+1):nrow(data),]  #scores for global climate
  scores.sp1 <- pca.cal$li[row.sp1,]					#scores for sp1
  scores.sp2 <- pca.cal$li[row.sp2,]					#scores for sp2
}
# calculation of occurence density and test of niche equivalency and similarity 
z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2,R=100)
}}

