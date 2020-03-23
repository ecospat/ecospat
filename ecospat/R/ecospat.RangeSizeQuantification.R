## This function estimates the range size of a species using standard IUCN criteria (AOO, EOO)
## or binary predictions of Species Distribution Models 

## FUNCTION'S ARGUMENTS
## bin.map:             Binary map (single layer or raster stack) from a Species Distribution Model.
## ocp:                 logical. Calculate occupied patch models from the binary map (predicted area occupied by the species).
## buffer:              numeric. Calculate occupied patch models from the binary map using a buffer (predicted area occupied by the species or within a buffer around the species).  
## eoo.around.model:    logical. The EOO around all positive predicted values from the binary map
## eoo.around.modelocp: logical. EOO around all positive predicted values of occupied patches.
## xy:                  xy-coordinates of the species presence
## EOO:                 logical. Extent of Occurrence. Convex Polygon around occurrences.
## Model.within.eoo:    logical. Postitive predicted Area within EOO.
## AOO:                 logical. Area of Occurrences.
## resol:               Resolution of the grid frame at which AOO should be calculated.
## AOO.circles:         logical. AOO calculated by circles around the occurrences instead of using a grid frame.
## d:                   Radius of the AOO.circles around the occurrences.
## lonlat:              Are these longitude/latidue coordinates? (Default FALSE)
## alphahull:           logical. Calculates the alpha-hull for the occurrences
## alpha:               numeric.
## return.obj:          logical. should the objects created to estimate range size be returned (rasterfiles and spatial polygons)
## save.obj:            logical. should opjects be saved on hard drive?
## save.rangesize:      logical. should range size estimations be saved on hard drive 
## directory:           directory toin which objects should be saved (Default: getwd())

## Details:

## Value
# 

## Author(s)
# Frank Breiner

## See Also
#   convHull, circles, ahull, ecospat.occupied.patch 


ecospat.rangesize <- function(bin.map = NULL,
                              ocp = TRUE,
                              buffer = 0,
                              eoo.around.model = TRUE,
                              eoo.around.modelocp = FALSE,
                              xy = NULL,
                              EOO = TRUE,
                              Model.within.eoo = TRUE,
                              AOO = TRUE,
                              resol = c(2000, 2000),
                              AOO.circles = FALSE,
                              d = sqrt((2000*2)/pi),
                              lonlat = FALSE,
                              # alpha.hull= F,
                              # alpha = 2,
                              return.obj = TRUE,
                              save.obj = FALSE,
                              save.rangesize = FALSE,
                              directory = getwd()){

  #if(ocp | EOO | alphahull | AOO & is.null(xy)){stop("need xy-coordinates of species occurrence")}
  #if(!requireNamespace(dismo)){stop("dismo package required!")}

  
  ### IUCN methods:
  if(Model.within.eoo | EOO){
    xy.eoo<-convHull(xy)
    eoo <- round(xy.eoo@polygons@polygons[[1]]@area) #quadrat km  
  }else{xy.eoo <- eoo <- NULL}  
  
  if(AOO){
    aoo.obj <- bin.map[[1]]
    res(aoo.obj) = resol
    aoo.obj[] <- 0
    
    aoo.obj <- rasterize(xy, aoo.obj, field=1) 
    aoo.rs <- as.numeric(freq( aoo.obj==1)[1,2]*prod(resol)) ### km2 AOO2
  }else{aoo.obj <- aoo.rs <- NULL}  
  
  if(AOO.circles){
#    circ   <- dismo::circles(xy,d=d,lonlat=lonlat)
    circ   <- circles(xy,d=d,lonlat=lonlat)
    circ.rs <- round(raster::area(circ@polygons)) 
  }else{circ.rs <- circ <- NULL}    
    
  # if(alpha.hull){
  #   
  #   #if(!requireNamespace(alphahull)){stop("alphahull package required!")}
  #   del<-delvor(xy)
  #   dv<-del$mesh
  #   mn <- mean(sqrt(abs(del$mesh[,3]-del$mesh[,5])^2+abs(del$mesh[,4]-del$mesh[,6])^2))*alpha
  #   h<-ahull(del,alpha=mn) 
  #   alpha.hull <- round(areaahull(h))
  # }else{h <- alpha.hull <- NULL}   

  
  ## Model range size estimation
  bin.map.rs <-NULL
  bin.map.ocp.rs <- bin.map.ocp <- NULL    
  eoo.around.mo <- list();eoo.around.mo.rs <- NULL    
  eoo.around.ocp <- NULL  
  eoo.around.mo.ocp <- list();eoo.around.mo.ocp.rs <- NULL    
  map.mo.within.eoo <- mo.within.eoo <- NULL
  
  if(!is.null(bin.map)){
  
    for(i in 1:nlayers(bin.map)){
      bin.m <- bin.map[[i]]
      bin.map.rs <- c(bin.map.rs,sum(bin.m[bin.m==1])*prod(res(bin.m)))
      names(bin.map.rs) <- names(bin.map[[1:i]])
      
      if(ocp | eoo.around.modelocp){
        if(i==1){
        bin.map.ocp  <- ecospat.occupied.patch(bin.map[[i]],xy, buffer=buffer)
        }else{
        bin.map.ocp  <- addLayer(bin.map.ocp, ecospat.occupied.patch(bin.map[[i]],xy, buffer=buffer))
        }
        bin.map.ocp.rs <- c(bin.map.ocp.rs,
        (sum(bin.map.ocp[[i]][bin.map.ocp[[i]]==2])/2)*prod(res(bin.map.ocp)))
        names(bin.map.ocp.rs) <- paste(names(bin.map[[1:i]]),"ocp",sep="_")      
      }

      if(eoo.around.model | eoo.around.modelocp){    
        ids <- init(bin.map[[i]], v='cell')
        # extract these
        cells <- extract(ids, ids[bin.map[[i]]==1])
        
        if(eoo.around.model){
        xy.model.eoo <- xyFromCell(ids,ids[bin.map[[i]]==1])
        eoo.around.mo[[i]] <- convHull(xy.model.eoo)
        eoo.around.mo.rs <- c(eoo.around.mo.rs,round(eoo.around.mo[[i]]@polygons@polygons[[1]]@area))
        names(eoo.around.mo.rs) <- paste(names(bin.map[[1:i]]),"eoo.around.model",sep="_")  
        names(eoo.around.mo) <- paste(names(bin.map[[1:i]]),"eoo.around.model",sep="_") 
        }       
      }
      if(eoo.around.modelocp){    
        xy.model.eoo.ocp <- xyFromCell(ids,ids[bin.map.ocp[[i]]==2])
        eoo.around.mo.ocp[[i]] <- convHull(xy.model.eoo.ocp)
        eoo.around.mo.ocp.rs <- c(eoo.around.mo.ocp.rs,round(eoo.around.mo.ocp[[i]]@polygons@polygons[[1]]@area))
        names(eoo.around.mo.ocp.rs) <- paste(names(bin.map[[1:i]]),"eoo.around.model.ocp",sep="_")  
        names(eoo.around.mo.ocp) <- paste(names(bin.map[[1:i]]),"eoo.around.model.ocp",sep="_") 
      }
        
      if(Model.within.eoo){
        d<-extract(bin.map[[i]], xy.eoo@polygons,cellnumbers = TRUE)
        mo.within.eoo <- c(mo.within.eoo,round(sum(d[[1]][,2],na.rm = TRUE)*prod(res(bin.map[[i]]))))
        
        cells <- d[[1]][,1][d[[1]][,2]==1 & !is.na(d[[1]][,2])]
        
        if(return.obj){
        if(i==1){
          map.mo.within.eoo <- bin.map[[i]]
          map.mo.within.eoo[cells] <-2
        }else{
          map.mo.within.eoo1 <- bin.map[[i]]
          map.mo.within.eoo1[cells] <-2
          map.mo.within.eoo  <- addLayer(map.mo.within.eoo, map.mo.within.eoo1)
        }
        }
        names(mo.within.eoo) <- paste(names(bin.map[[1:i]]),"mo.within.eoo",sep="_")
      }

    }
  }
  
  range.sizes <-list(AOO = aoo.rs, 
  AOO.circle = circ.rs, 
  # alpha.hull = alpha.hull, 
  EOO = eoo, model=bin.map.rs, 
  models.ocp= bin.map.ocp.rs, 
  eoo.around.model = eoo.around.mo.rs, 
  eoo.around.mo.ocp =eoo.around.mo.ocp.rs, 
  model.within.eoo = mo.within.eoo)
  objects <- list(AOO = aoo.obj, 
  AOO.circle = circ, 
  # alpha.hull = h, 
  EOO = xy.eoo, 
  models.ocp = bin.map.ocp,
  eoo.around.model = eoo.around.mo, 
  eoo.around.mo.ocp = eoo.around.mo.ocp, 
#  =model.within.eoo = map.mo.within.eoo)
  model.within.eoo = map.mo.within.eoo)
  

  if(return.obj){
    return(list(RangeSize= range.sizes, RangeObjects=objects))
  }else{
    return(list(RangeSize= range.sizes))   
  }
  
  
  if(save.obj){
    save(objects, file=paste(directory,"rangesize_objects",sep="/"))
  }
  if(save.obj){
    save(save.rangesize, file=paste(directory,"rangesize",sep="/"))
  }
}


### get occupied patches

## This function estimates the range.size of a species using standard IUCN criteria (AOO, EOO)
## or using Species Distribution Models 

## FUNCTION'S ARGUMENTS
## bin.map:             Binary map (single layer or raster stack) from a Species Distribution Model.
## Sp.occ.xy:           xy-coordinates of the species presence
## buffer:              numeric. Calculate occupied patch models from the binary map using a buffer (predicted area occupied by the species or within a buffer around the species, for details see ?extract).  

## Details:

## Value
# a RasterLayer with value 1 ()

## Author(s)
# Frank Breiner



ecospat.occupied.patch <- function(bin.map, Sp.occ.xy, buffer = 0){
  cl <- clump(bin.map)
  d <- raster::extract(cl,Sp.occ.xy,buffer=buffer)
  d <- unique(na.omit(unlist(d)))
  b.map <- raster::match(cl,d) + bin.map
  names(b.map) <- names(bin.map)
  return(b.map)}

