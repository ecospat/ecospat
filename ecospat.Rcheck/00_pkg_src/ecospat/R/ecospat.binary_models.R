### Author: Ruben G Mateo
### Pred}{Predicted suitability values (from 0 to 1000). A RasterStack object containing models predictions ( Output from biomod2 in raster format).}
### Sp.occ.xy}{Ocurrences of the species. A dataframe object with two columns: longitude and latitude. Coordinate systems other than longitude and latitude can be used, for example "x" and "y".}
### Percentage}{The percentage of omission error used to generate the binary model.}
#ecospat.binary.model <-function(Pred, Sp.occ.xy, Percentage)
#{
#  # Calculate threshold
#  #####################
#  Pres <- extract (Pred, Sp.occ.xy)
#  Pres.sort <- data.frame(sort(Pres))
#  Num.pres <- nrow(Sp.occ.xy)
#  Num <- as.integer(((Num.pres * Percentage)/100) + 0.5)
#  #ifelse (Num == 1, Num <- 2, Num)
#  Threshold <- sort(Pres) [Num]
#
#
#  # Generation of binary model
#  ############################
#  Threshold.Com <- c(0, Threshold, 0, Threshold, 1000, 1)
#  Threshold.Com.b <- matrix (Threshold.Com, ncol=3, byrow=TRUE)
#  Pred.binary <- reclassify (Pred,Threshold.Com.b)
#}


## Pred:{Predicted suitability values (from 0 to 1000). A RasterStack object containing models predictions ( Output from biomod2 in raster format).}
## Threshold:{}


### New version with contributions of Frank Breiner

ecospat.binary.model <- function (Pred, Threshold)
{
  Threshold.Com <- c(0, Threshold, 0, Threshold, maxValue(Pred), 1)
  Threshold.Com.b <- matrix(Threshold.Com, ncol = 3, byrow = T)
  Pred.binary <- reclassify(Pred, Threshold.Com.b)
  names(Pred.binary) <- names(Pred)
  return(Pred.binary)
}