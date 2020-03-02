
## Written by Olivier Broennimann. Departement of Ecology and Evolution (DEE) & Institude of Earth Surface Dynamics
## University of Lausanne. Switzerland. April 2018.
##
## DESCRIPTION
##
## functions to project niche quantification (objets z calculated with ecospat.grid.clim.dyn) onto the geographical space
## 
## list of functions:
##
## ecospat.niche.zProjGeo(z1,env,cor)
## projects the density of occurrence in space 
## z1 is an object created by ecospat.grid.clim.dyn
## env is a RasterStack of environmental variables corresponding to the background (glob in ecospat.grid.clim.dyn)
## cor tells if the corrected or uncorrected occurrence density should be projected

## ecospat.niche.zProjGeo(z1,z2,env,index)
## projects the dynamic indexes ("stability", "unfilling" and "expansion" in space
## z1 and z2 are objects created by ecospat.grid.clim.dyn
## env is a RasterStack of environmental variables corresponding to the background (glob in ecospat.grid.clim.dyn)
## index tells which which index to project ("stability", "unfilling" or "expansion")

ecospat.niche.zProjGeo <- function(z1,env,cor=FALSE){

  XY <- rasterToPoints(env)[,1:2] #geographical coordinates of each point of the background

  if (cor==FALSE) Z1<-extract(z1$z.uncor,z1$glob1) # occurrence density (niche) for each point of the background
  if (cor==TRUE) Z1<-extract(z1$z.cor,z1$glob1)
  XYZ1<-cbind(XY,Z1)
  geoz1<-rasterFromXYZ(XYZ1)

  return(geoz1)
}

#########################

ecospat.niche.dynIndexProjGeo <- function(z1,z2,env,index=NULL){

  XY <- rasterToPoints(env)[,1:2] #geographical coordinates of each point of the background

  if(index=="stability") ind<-raster(ecospat.niche.dyn.index(z1,z2)$dyn==2)
  if(index=="unfilling") ind<-raster(ecospat.niche.dyn.index(z1,z2)$dyn==-1) 
  if(index=="expansion") ind<-raster(ecospat.niche.dyn.index(z1,z2)$dyn==1) 
  if(index!="stability"&index!="unfilling"&index!="expansion") stop("set index as stability,unfilling or expansion")
  ind@extent <- z1$z.uncor@extent

  ZS<-extract(ind,z1$glob)
  XYZS<-cbind(XY,ZS)
  geozS<-rasterFromXYZ(XYZS)

  return(geozS)
}

