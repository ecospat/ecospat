\name{ecospat.makeDataFrame}

\alias{ecospat.makeDataFrame}

\title{Make Data Frame}

\description{Create a biomod2-compatible dataframe. The function also enables to remove duplicate presences within a pixel and to set a minimum distance between presence points to avoid autocorrelation. Data from GBIF can be added.}

\usage{ecospat.makeDataFrame (spec.list, expl.var, use.gbif=FALSE, precision=NULL,
year=NULL, remdups=TRUE, mindist=NULL, n=1000, type='random', PApoint=NULL,
ext=expl.var, tryf=5)}

\arguments{
  \item{spec.list}{Data.frame or Character. The species occurrence information must be a data.frame in the form:
  \'x-coordinates\' , \'y-coordinates\' and \'species name\' (in the same projection/coordinate system as expl.var!).}
  \item{expl.var}{a RasterStack object of the environmental layers.}
  \item{use.gbif}{Logical. If TRUE presence data from GBIF will be added. It is also possible to use GBIF data only.
   Default: FALSE. See ?gbif {dismo} for more information. Settings: geo=TRUE, removeZeros=TRUE,
   all sub-taxa will be used.
   \'species name\' in spec.list must be in the form: \'genus species\',
   \'genus_species\' or \'genus.species\'. If there is no species information available on GBIF an error is returned.
   Try to change species name (maybe there is a synonym) or switch use.gbif off.}
  \item{precision}{Numeric. Use precision if use.gbif = TRUE to set a minimum precision of
   the presences which should be added. For precision = 1000
   e.g. only presences with precision of at least 1000 meter will be added from GBIF.
   When precision = NULL all presences from GBIF will be used, also presences
   where precision information is NA.}
  \item{year}{Numeric. Latest year of the collected gbif occurrences. If year=1960 only occurrences which were collected since 1960 are used.}
  \item{remdups}{Logical. If TRUE (Default) duplicated presences within a raster pixel
   will be removed. You will get only one presence per pixel.}
  \item{mindist}{Numeric. You can set a minimum distance between presence points to avoid
   autocorrelation. nndist {spatstat} is used to calculate the nearest neighbour (nn)
   for each point. From the pair of the minimum nn, the point is removed of which the second
   neighbour is closer. Unit is the same as expl.var.}
  \item{n}{number of Pseudo-Absences. Default 1000.}
  \item{type}{sampling dessign for selecting Pseudo-Absences. If \'random\' (default) background points are selected with the 
   function randomPoints {dismo}. When selecting another sampling type (\'regular\', \'stratified\', \'nonaligned\',
   \'hexagonal\', \'clustered\' or \'Fibonacci\') spsample {sp} will be used. This can immensely increase computation time and RAM usage if ext is a raster,
   especially for big raster layers because it must be converted into a  \'SpatialPolygonsDataFrame\' first.}
  \item{PApoint}{data.frame or SpatialPoints. You can use your own set of Pseudo-Absences
   instead of generating new PAs. Two columns with \'x\' and \'y\' in the same
   projection/coordinate system as expl.var!}
  \item{ext}{a Spatial Object or Raster object. Extent from which PAs should be selected from (Default is expl.var).} 
  \item{tryf}{numeric > 1. Number of trials to create the requested Pseudo-Absences after removing NA points (if type='random'). See ?randomPoints {dismo}}
}

\details{If you use a raster stack as explanatory variable and you want to model many species in a loop with Biomod, formating data will last very long as presences and PA's have to be extracted over and over again from the raster stack. To save computation time, it is better to convert the presences and PAs to a data.frame first.}

\value{A data.frame object which can be used for modeling with the Biomod package.}

\author{Frank Breiner \email{frank.breiner@unil.ch}}

\examples{
\dontrun{
files <- list.files(path=paste(system.file(package="dismo"),
                               '/ex', sep=''), pattern='grd', full.names=TRUE )
predictors <- raster::stack(files[c(9,1:8)])   #file 9 has more NA values than
# the other files, this is why we choose it as the first layer (see ?randomPoints)

file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file, header=TRUE, sep=',')[,c(2,3,1)]
head(bradypus)

random.spec <- cbind(as.data.frame(randomPoints(predictors,50,extf=1)),species="randomSpec")
colnames(random.spec)[1:2] <- c("lon","lat")

spec.list <- rbind(bradypus, random.spec)

df <- ecospat.makeDataFrame(spec.list, expl.var=predictors, n=5000)
head(df)

plot(predictors[[1]])
points(df[df$Bradypus.variegatus==1, c('x','y')])
points(df[df$randomSpec==1, c('x','y')], col="red")

}
}
