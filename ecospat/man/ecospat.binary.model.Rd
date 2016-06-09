\name{ecospat.binary.model}

\alias{ecospat.binary.model}

\title{Generate Binary Models}

\description{Generate a binary map from a continuous model prediction (i.e., values from 0 to 1000).}

\usage{ecospat.binary.model (Pred, Sp.occ.xy, Percentage)}

\arguments{
  \item{Pred}{Predicted suitability values (from 0 to 1000). A RasterStack object containing models predictions ( Output from biomod2 in raster format).}
  \item{Sp.occ.xy}{Ocurrences of the species. A dataframe object with two columns: longitude and latitude. Coordinate systems other than longitude and latitude can be used, for example "x" and "y".}
  \item{Percentage}{The percentage of omission error used to generate the binary model.}
}

\details{This function generates a binary model prediction (presence/absence) from an original model (continuous values from 0 to 1000) applying a threshold of maximum acceptable error of false negatives (i.e. percentage of the presence predicted as absences, omission error).}

\value{The binary model prediction (presence/absence) generated.}

\references{Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation}, \bold{24}: 38-49.

Engler, R., A Guisan and L. Rechsteiner. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. \emph{Journal of Applied Ecology}, \bold{41}, 263-274.

Guisan, A., O. Broennimann, R. Engler, M. Vust, N.G. Yoccoz, A. Lehmann and N.E. Zimmermann. 2006. Using niche-based models to improve the sampling of rare species. \emph{Conservation Biology}, \bold{20}, 501-511.

}

\author{Ruben G. Mateo \email{rubeng.mateo@gmail.com}}

\examples{

\dontrun{

#Run biomod2 to produce a model prediction 
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv", package="biomod2"))

myRespName <- 'GuloGulo'
# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])
# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
# load the environmental raster layers (could be .img, ArcGIS
# rasters or any supported format by the raster package)
# Environmental variables extracted from Worldclim (bio_3, bio_4,
# bio_7, bio_11 & bio_12)
myExpl = stack( system.file( "external/bioclim/current/bio3.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))


myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
myBiomodData
myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM'),
  models.options = myBiomodOption,
  NbRunEval=1,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

myBiomodModelOut
myBiomodModelEval <- get_evaluations(myBiomodModelOut)

myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = TRUE,
  prob.cv = TRUE,
  prob.ci = TRUE,
  prob.ci.alpha = 0.05,
  prob.median = TRUE,
  committee.averaging = TRUE,
  prob.mean.weight = TRUE,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM

myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = FALSE,
  output.format = '.grd')

myBiomodProj
plot(myBiomodProj, str.grep = 'GLM')
#
Pred <-get_predictions(myBiomodProj)
Sp.occ.xy <- DataSpecies[DataSpecies[,5]==1,2:3]
Percentage <- 7

binary.model<-ecospat.binary.model (Pred, Sp.occ.xy, Percentage)
plot(binary.model)
}
}