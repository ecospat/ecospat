\name{ecospat.rcls.grd}
\alias{ecospat.rcls.grd}

\title{Reclassifying grids function}

\description{Function for reclassifying grid files to get a combined statification from more than one grid}

\usage{ecospat.rcls.grd(in_grid,no.classes)
}
\arguments{
    \item{in_grid}{The grid to be reclassified.}
    \item{no.classes}{The number of desired new classes.}
}
\details{
This function reclassifies the input grid into a number of new classes that the user defines. The boundaries of each class are decided automatically by splitting the range of values of the input grid into the user defined number of classes.
}

\value{
Returns a reclassified Raster object
}
\author{
Achilleas Psomas \email{achilleas.psomas@wsl.ch} and Niklaus E. Zimmermann \email{niklaus.zimmermann@wsl.ch}
}

\examples{

   \dontrun{
bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))
B3.rcl<-ecospat.rcls.grd(bio3,9) 
B12.rcl<-ecospat.rcls.grd(bio12,9)
B3B12.comb <- B12.rcl+B3.rcl*10

# Plotting a histogram of the classes
hist(B3B12.comb,breaks=100,col=heat.colors(88)) 
# Plotting the new RasterLayer (9x9 classes)
plot(B3B12.comb,col=rev(rainbow(88)),main="Stratified map") 

}}
