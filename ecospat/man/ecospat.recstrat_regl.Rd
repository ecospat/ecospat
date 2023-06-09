\name{ecospat.recstrat_regl}
\alias{ecospat.recstrat_regl}

\title{
  Random Ecologically Stratified Sampling of equal numbers
}
\description{
  This function randomly takes an equal number of samples per class in the stratification layer. 
}
\usage{
  ecospat.recstrat_regl(in_grid, sample_no)
}
\arguments{
  \item{in_grid}{
    The stratification grid (SpatRaster) to be sampled.
  }
  \item{sample_no}{
    The total number of pixels to be sampled.
  }
}
\details{
  The number of classes in the stratification layer is determined automatically from the integer input map. 
  If the number of pixels in a class is higher than the number of samples, 
  then a random selection without re-substitution is performed, 
  otherwise all pixels of that class are selected.
}
\value{
  Returns a dataframe with the selected sampling locations their coordinates and the strata they belong in.
}
\author{
  Achilleas Psomas \email{achilleas.psomas@wsl.ch} and Niklaus E. Zimmermann \email{niklaus.zimmermann@wsl.ch}
}

\seealso{
  \code{\link{ecospat.recstrat_prop}}
  \code{\link{ecospat.rcls.grd}}
}
\examples{

  library(terra)
  library(classInt)
  library(biomod2)
  
  data("bioclim_current")
  bioclim_current <- rast(bioclim_current)
  bio3 <- bioclim_current[["bio3"]]
  bio12 <- bioclim_current[["bio12"]]

    
  B3.rcl<-ecospat.rcls.grd(bio3,9) 
  B12.rcl<-ecospat.rcls.grd(bio12,9)
  B3B12.comb <- B12.rcl+B3.rcl*10
    
  B3B12.regl_samples <- ecospat.recstrat_prop(B3B12.comb,100)
  
  plot(B3B12.comb)
  points(B3B12.regl_samples$x,B3B12.regl_samples$y,pch=16,cex=0.6,col=B3B12.regl_samples$class)
}
