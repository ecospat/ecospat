
## Written by Olivier Broennimann. Departement of Ecology and Evolution
## (DEE) & Institude of Earth Surface Dynamics
## University of Lausanne. Switzerland. April 2018.
##
## ecospat.niche.dynIndexProjGeo updated by Tyler Smith, April 2023
##
## DESCRIPTION
##
## functions to project niche quantification (objets z calculated with
## ecospat.grid.clim.dyn) onto the geographical space
## 
## list of functions:
##
## ecospat.niche.zProjGeo(z1,env,cor)
## projects the density of occurrence in space 
## z1 is an object created by ecospat.grid.clim.dyn
## env is a RasterStack of environmental variables corresponding to the background (glob in ecospat.grid.clim.dyn)
## cor tells if the corrected or uncorrected occurrence density should be projected

## ecospat.niche.dynIndexProjGeo(z1, z2, env)
## projects the dynamic indexes ("stability", "unfilling" and "expansion")
## in space
## z1 and z2 are objects created by ecospat.grid.clim.dyn
## env is a RasterStack of environmental variables corresponding to the
## background (glob in ecospat.grid.clim.dyn)

ecospat.niche.zProjGeo <- function(z1,env,cor=FALSE){

  XY <- rasterToPoints(env)[,1:2] #geographical coordinates of each point of the background

  if (cor==FALSE) Z1<-extract(z1$z.uncor,z1$glob1) # occurrence density (niche) for each point of the background
  if (cor==TRUE) Z1<-extract(z1$z.cor,z1$glob1)
  XYZ1<-cbind(XY,Z1)
  geoz1<-rasterFromXYZ(XYZ1)

  return(geoz1)
}

#########################

ecospat.niche.dynIndexProjGeo <- function(z1, z2, env) {
  zRast <- z1$w + 2*z2$w

  ## We need **all** the cells from the env raster, even the NA ones. Those
  ## are not present in z1$glob or z2$glob, so we need the original env
  ## raster.

  ## all env cells:
  envVals <- raster::getValues(env)[,1]

  ## index values for the non-NA cells:
  zCats <- raster::extract(zRast, z1$glob)

  ## add index values to the non-NA cells:
  envVals[!is.na(envVals)] <- zCats
  
  ## create raster (without using values() which triggers namespace conflicts)

  zCatRast<-raster(matrix(envVals, nrow=env[[1]]@nrows, ncol=env[[1]]@ncols,byrow = TRUE),
                   xmn=env@extent[1],xmx=env@extent[2],ymn=env@extent[3],ymx=env@extent[4])
  
  return(zCatRast)
}
