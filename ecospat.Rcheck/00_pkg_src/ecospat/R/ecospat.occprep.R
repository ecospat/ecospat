## Remove species occurrences in a dataframe that are closer to each other than a specified distance threshold.

## FUNCTION'S ARGUMENTS
## xy:        A dataframe with xy-coordinates (x-column must be named 'x' and y-column 'y')
## min.dist:  The minimun distance between points in the sub-dataframe
## by:        Grouping element in the dataframe (e.g. species)

## Details:

## Value
# A subset of df with the columns specified in colvar.

## Author(s)
# Frank Breiner frank.breiner@unil.ch with contributions of Olivier Broennimann


ecospat.occ.desaggregation <- function(xy, min.dist,by=NULL){
  if(is.null(xy$x)|is.null(xy$y)){
    stop("no x and/or y column")
  }

  if(!is.null(by)){
    xy[,which(names(xy)==by)] <- factor(xy[,which(names(xy)==by)])
  }
  new.data <- NULL

  if(is.null(by)){
    to<-1
  }else{
    to<-nlevels(xy[,which(names(xy)==by)])
  }
  for(i in 1:to){
    if(to>1){
    print(paste("desaggregate species",i))}
    if(!is.null(by)){
      del.min.dist <- xy[xy[,by]==levels(xy[,by])[i],]
    }else{
        del.min.dist <- xy
    }

    if(sum(duplicated(paste(del.min.dist$x,del.min.dist$y)))>0){
    stop(paste("duplicated values",levels(xy[,by])[i],sep=" "))
    }

    repeat{
      
    nn1 <- nndist(del.min.dist[,"x"],del.min.dist[,"y"])       # calculate distance nearest neighbour
    if (sum(nn1 < min.dist) == 0){
      break
    }
    # iteratively removing points starting with the one having the minimal distance to the nearest neighbour
    nn2 <- nndist(del.min.dist[,"x"],del.min.dist[,"y"],k=2)

        del1 <- nn1 == min(nn1) 
        del2 <- nn2==min(nn2[del1])
        delk <- del1 & del2 
        if(sum(del2)>1){
          for(k in 3:8){

          nn <- nndist(del.min.dist[,"x"],del.min.dist[,"y"],k=k)
          delk <- delk & nn==min(nn[delk])
          if(sum(nn[delk]==min(nn[delk]))>1){break}
          }
        }
            
        # from the two points which are the nearest neighbours of the whole set, remove the one closest to the second neighbour
        del.min.dist <- del.min.dist[-(which(delk)[1]),]
    }

    new.data <- rbind(new.data,del.min.dist)
  }
  result <- list(initial = nrow(xy), kept = nrow(new.data), out = nrow(xy)-nrow(new.data))
  print(result)
  return(xy=new.data[,!colnames(new.data) %in% c("nn","nn2","id")])
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##
## add environmental values to a species dataframe.
## the xy (lat/long) coordinates of the species occurrences are compared to those of the environment dataframe
## and the value of the closest pixel is added to the species dataframe. 
## when the closest environment pixel is more distant than resolution, NA is added instead of the value.
## (similar to sample() in ArcGIS)

##ARGUMENTS
##dfsp: species dataframe with x, y and optional other variables
##colspxy: the range of columns for x and y in dfsp
##colspkept: the columns of dfsp that should be kept in the final dataframe (by default: xy )
##dfvar: environmental dataframe with x, y and environmental variables
##colvarxy: the range of columns for x and y in dfvar
##colvar: the range of enviromental variables columns in dfvar. (by default: all exept xy )
##resolution: distance between x,y of species and environmental datafreme after which values shouldn't be added 
##(typically, the resolution of the data in dfvar)

ecospat.sample.envar <- function(dfsp, colspxy, colspkept = "xy", dfvar, colvarxy, colvar = "all",
  resolution) {

  if (sum(colspkept == "xy") == 1)
    colspkept <- colspxy
  if (sum(colvar == "all") == 1) {
    if (!is.null(colspkept))
      colvar <- (1:ncol(dfvar))[-colvarxy]
    if (is.null(colspkept))
      colvar <- (1:ncol(dfvar))
  }
  colspx <- colspxy[1]
  colspy <- colspxy[2]
  colvarx <- colvarxy[1]
  colvary <- colvarxy[2]

  x <- dfsp[, colspx]
  X <- dfvar[, colvarx]
  y <- dfsp[, colspy]
  Y <- dfvar[, colvary]

  train <- data.frame(matrix(nrow = nrow(dfsp), ncol = length(colvar)))
  names(train) <- names(dfvar)[colvar]

  for (i in 1:nrow(dfsp)) {
    dist <- sqrt((X - x[i])^2 + (Y - y[i])^2)
    min <- min(dist)
    if (min <= resolution) {
      if (length(colvar) > 1)
        train[i, ] <- dfvar[dist == min, colvar][1, ]
      if (length(colvar) == 1)
        train[i, ] <- dfvar[dist == min, colvar][1]
    }
  }


  if (!is.null(colspkept))
    final <- cbind(dfsp[, colspkept], train)
  if (is.null(colspkept))
    final <- train

  return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##Investigate spatial autocorrelation by drawing a mantel Correlogram (autocorrelation vs distance)
##
##ARGUMENTS
##df: dataframe with x, y, and variables
##colxy: the range of columns for x and y in df
##colvar: the range of columns for variables in df
##n: number of random occurences used for the test (computation time increase tremendiously when using more than 500occ.) 	
##max: maximum distance to be computed in the correlogram
##nclass: number of class of distance to be computed in the correlogram
##nperm: number of permutation in the randomization process



ecospat.mantel.correlogram <- function(dfvar, colxy, n, colvar, max, nclass, nperm) {
  envnorm <- data.frame(t((t(dfvar[, colvar]) - apply(dfvar[, colvar], 2, mean))/apply(dfvar[, colvar],
    2, sd)))
  row.rand <- sample(1:nrow(dfvar), n, replace = TRUE)
  envdist <- dist(envnorm[row.rand, ])
  geodist <- dist(dfvar[row.rand, colxy])
  b <- seq(from = min(geodist), to = max, length.out = nclass)
  crlg <- mgram(envdist, geodist, breaks = b, nperm = nperm)
  plot(crlg)
  abline(h = 0)
}



##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##randomly sample pseudoabsences from an environmental dataframe covering the study area
##A minimum distance from presences can be set.
##ARGUMENTS
##nbabsences: number of pseudoabsences desired 
##glob: environmental dataframe covering the study area to sample, with x,y 
##colxyglob: the range of columns for x and y in glob
##colvar: the range of columns for x and y in glob. colvar="all" keeps all the variables in glob in the final dataframe. colvar=NULL keeps only x and y
##presence: occurence dataframe 
##colxypresence: the range of columns for x and y in presence
##mindist: minimum distance from prensences closer to wich pseudoabsences shouldn't be drawn (buffer distance around presences)

ecospat.rand.pseudoabsences <- function(nbabsences, glob, colxyglob, colvar = "all", presence, colxypresence,
  mindist) {

  colxglob <- colxyglob[1]
  colyglob <- colxyglob[2]
  colxpresence <- colxypresence[1]
  colypresence <- colxypresence[2]

  keep <- c()

  no.i <- 1
  while (no.i <= nbabsences) {
    ki <- sample(1:nrow(glob), 1)
    if (sum(((glob[ki, colxglob] - presence[, colxpresence])^2 + (glob[ki, colyglob] - presence[,
      colypresence])^2) <= mindist^2) == 0) {
      keep[no.i] <- ki
      no.i <- no.i + 1
    }
  }
  if (sum(colvar == "all") == 1)
    colvar <- (1:ncol(glob))[-colxyglob]
  if (!is.null(colvar))
    pseudoabs <- glob[keep, c(colxyglob, colvar)]
  if (is.null(colvar))
    pseudoabs <- glob[keep, colxyglob]

  return(pseudoabs)
}