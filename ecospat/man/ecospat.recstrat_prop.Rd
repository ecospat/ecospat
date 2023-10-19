\name{ecospat.recstrat_prop}
\alias{ecospat.recstrat_prop}

\title{
  Random Ecologically Stratified Sampling of propotional numbers
}
\description{
  This function randomly collects a user-defined total number of samples from the stratification layer. 
}
\usage{
  ecospat.recstrat_prop(in_grid, sample_no)
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
  The number of samples per class are determined proportional to the abundance of each class.
  The number of classes in the stratification layer are determined automatically from the integer input map.
  If the proportion of samples for a certain class is below one then no samples are collected for this class.
}
\value{
  Returns a dataframe with the selected sampling locations their coordinates and the strata they belong in.
}
\author{
  Achilleas Psomas \email{achilleas.psomas@wsl.ch} and Niklaus E. Zimmermann \email{niklaus.zimmermann@wsl.ch}
}

\seealso{
  \code{\link{ecospat.recstrat_regl}}
  \code{\link{ecospat.rcls.grd}}
}
\examples{
 
    library(terra)
    library(classInt)
    library(biomod2)
    
    data("bioclim_current")
    bioclim_current <- terra::rast(bioclim_current)
    bio3 <- bioclim_current[["bio3"]]
    bio12 <- bioclim_current[["bio12"]]

    
    B3.rcl<-ecospat.rcls.grd(bio3,9) 
    B12.rcl<-ecospat.rcls.grd(bio12,9)
    B3B12.comb <- B12.rcl+B3.rcl*10
    
    B3B12.prop_samples <- ecospat.recstrat_prop(B3B12.comb,100)
    
    plot(B3B12.comb)
    points(B3B12.prop_samples$x,B3B12.prop_samples$y,pch=16,cex=0.6,col=B3B12.prop_samples$class)
}
