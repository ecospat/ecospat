
## Written by Olivier Broennimann. Departement of Ecology and Evolution
## (DEE) & Institude of Earth Surface Dynamics
## University of Lausanne. Switzerland. April 2018.
## 
## ecospat.niche.dynIndexProjGeo updated by Tyler Smith, April 2023 and Flavien Collart
##
## DESCRIPTION
##
## functions to project niche quantification (objets z calculated with
## ecospat.grid.clim.dyn) onto the geographical space
## 
## list of functions:
##
## ecospat.niche.zProjGeo(z,zproj=NULL,env,cor=FALSE)
## projects the density of occurrence in space 
## z and zproj are objects created by ecospat.grid.clim.dyn
## env is a RasterStack of environmental variables corresponding to the
## background of z if zproj=NULL or zproj if not NULL (glob1 in ecospat.grid.clim.dyn)
## cor tells if the corrected or uncorrected occurrence density should be projected
##
## ecospat.niche.dynIndexProjGeo(z1,z2,proj=0,env)
## projects the dynamic indexes ("stability", "unfilling" and "expansion") in space
## z1 and z2 are objects created by ecospat.grid.clim.dyn
## proj indicate the type of projection: 0 for the global range, 1 for the range of z1, 2 for the range of z2
## env is a SpatRaster environmental variables corresponding to the background (glob in ecospat.grid.clim.dyn
## if proj = 0, glob1 if proj = 1 or proj = 2)

ecospat.niche.zProjGeo <- function(z,zproj=NULL,env,cor=FALSE){
  XY <- terra::crds(terra::as.points(env)) #geographical coordinates of each point of the background
  
  if(is.null(zproj)){
    if(cor==FALSE) Z<-extract(z$z.uncor,z$glob1) # occurrence density (niche) for each point of the background
    if(cor==TRUE) Z<-extract(z$z.cor,z$glob1)
  }else{
    if(cor==FALSE) Z<-extract(z$z.uncor,zproj$glob1) # occurrence density (niche) for each point of the background
    if(cor==TRUE) Z<-extract(z$z.cor,zproj$glob1)
  }
  XYZ<-cbind(XY,Z)
  geoz<-terra::rast(XYZ, type="xyz")
  
  return(geoz)
}

#########################

ecospat.niche.dynIndexProjGeo <- function(z1, z2, proj=0, env) {
  zRast <- z1$w + 2*z2$w
  #0: outside of both niches
  #1: only in z1 niche-> unfilling
  #2: only in z2 niche-> expansion
  #3: in both z1 and z2 niches -> stability

  ## We need **all** the cells from the env raster, even the NA ones. Those
  ## are not present in z1$glob or z2$glob, so we need the original env
  ## raster.

  ## all env cells:
  envVals <- terra::values(env)[,1]

  ## index values for the non-NA cells:
  if(proj==0){ zCats <- terra::extract(zRast, z1$glob)[,1]}
  if(proj==1){ zCats <- terra::extract(zRast, z1$glob1)[,1]}
  if(proj==2){ zCats <- terra::extract(zRast, z2$glob1)[,1]}

  ## add index values to the non-NA cells:
  envVals[!is.na(envVals)] <- zCats
  
  ## create grid (without using values() which triggers namespace conflicts)

  zCatRast<-terra::rast(matrix(envVals, nrow=nrow(env[[1]]), ncol=ncol(env[[1]]),byrow = TRUE),
                        extent = terra::ext(env))
  
  return(zCatRast)
}
