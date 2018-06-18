# @ ecospat.rcls.grd This function was written by Achilleas Psomas for reclassifying grids
# The purpose is to get a combined statification from more than one grid
# written by achilleas.psomas@wsl.ch

### THIS FUNCTION REQUIRES THAT THE PACKAGE 'classInt' is installed

ecospat.rcls.grd <- function(in_grid, no.classes) {
  new_classes <- classIntervals(getValues(in_grid), no.classes, style = "equal")
  # new_classes <-classIntervals(getValues(in_grid), no.classes , style = 'sd')
  new_classes_breaks <- new_classes$brks
  new_classes_limits <- matrix(ncol = 3, nrow = no.classes)
  classes <- 1:no.classes
  for (i in 1:no.classes) {
    new_classes_limits[i, 1] <- new_classes_breaks[i]
    new_classes_limits[i, 2] <- new_classes_breaks[i + 1]
    new_classes_limits[i, 3] <- classes[i]
  }
  in_grid_reclass <- reclassify(in_grid, new_classes_limits, include.lowest = TRUE)
  return(in_grid_reclass)
}


# ---------------------------------------------------------------------------#
# @ ecospat.recstrat_prop                                                    #
# Purpose: Random Ecologically Stratified Sampling of propotional numbers    #
# Authors: Achilleas Psomas, Niklaus E. Zimmermann                           #
# Date:    23.9.2009                                                   	     #
# Contact: achilleas.psomas(at)wsl.ch, nez(at)wsl.ch                   	     #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  	     #
# This function randomly collects a user-defined total number of samples     #
# from the stratification layer. The number of samples per class             #
# are determined proportional to the abundance of each class.                #
# The number of classes in the stratification layer are determined 	         #
# automatically from the integer input map.                                  #
# If the proportion of samples for a certain class is below 1 then no        #
# samples are collected for this class.                                      #
# REQUIRED INPUT:                                                      	     #
#   in_grid   = The stratification grid to be sampled                        #
#   sample_no = The total number of pixels to be sampled                     #
# ---------------------------------------------------------------------------#

ecospat.recstrat_prop <- function(in_grid, sample_no) {
  strata <- na.omit(unique(getValues(in_grid)))
  strata_no <- length(strata)
  in_grid_SPixels <- as(in_grid, "SpatialPointsDataFrame")
  total_pixels <- nrow(in_grid_SPixels@data)

  strata_stats <- table(in_grid_SPixels$layer)

  strata_stats_sorted <- as.data.frame(sort(strata_stats, decreasing = TRUE))
  pixels_largest_strata <- max(strata_stats_sorted$Freq)
  proportion_largest_strata <- round((pixels_largest_strata * sample_no)/total_pixels)

  result_list <- list()
  for (j in 1:length(strata)) {
    grid_sel <- in_grid_SPixels[in_grid_SPixels@data[, 1] == strata[j], ]
    proportion <- ceiling(log(dim(grid_sel@data)[1])/log(pixels_largest_strata) * proportion_largest_strata)

    optimal_samples_per_class <- ifelse(proportion < dim(grid_sel@data)[1], proportion, dim(grid_sel@data)[1])
    sp_points <- grid_sel[sample(1:nrow(grid_sel), optimal_samples_per_class, replace = FALSE), ]
    sample_points <- cbind(sp_points@coords, class = sp_points@data[, 1])
    result_list[[j]] <- sample_points
  }
  result <- data.frame(do.call("rbind", result_list))
  names(result) <- c("x", "y", "class")
  return(result)
}


  # -------------------------------------------------------------------- #
  # @ recstrat_equal                                                     #
  # Purpose: Random Ecologically Stratified Sampling of equal numbers    #
  # Authors: Achilleas Psomas, Niklaus E. Zimmermann                     #
  # Date:    23.9.2009                                                   #
  # Contact: achilleas.psomas(at)wsl.ch, nez(at)wsl.ch                   #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
  # This function randomly takes an equal number of samples per class    #
  # in the stratification layer. The number of classes in the stratifi-  #
  # cation layer is determined automatically from the integer input map. #
  # If the number of pixels in a class is higher than the number of sam- #
  # ples, then  a random selection without re-substitution is performed, #
  # otherwise all pixels of that class are selected.                     #
  # REQUIRED INPUT:                                                      #
  #   in_grid   = The stratification grid to be sampled                  #
  #   sample_no = The total number of pixels to be sampled               #
  # -------------------------------------------------------------------- #

ecospat.recstrat_regl <- function(in_grid, sample_no) {
  strata <- na.omit(unique(getValues(in_grid)))
  strata_no <- length(strata)
  if (sample_no < strata_no) {
    stop("Stoping Execution: The number of samples is lower the total number of unique classes")
  } else {
    samples_per_class <- round(sample_no/strata_no)
  }

  in_grid_SPixels <- as(in_grid, "SpatialPointsDataFrame")

  result_list <- list()
  for (j in 1:length(strata)) {
    grid_sel <- in_grid_SPixels[in_grid_SPixels@data[, 1] == strata[j], ]
    optimal_samples_per_class <- ifelse(dim(grid_sel@data)[1] > samples_per_class, samples_per_class,
      dim(grid_sel@data)[1])
    sp_points <- grid_sel[sample(1:nrow(grid_sel), optimal_samples_per_class, replace = FALSE), ]
    result_list[[j]] <- cbind(sp_points@coords, sp_points@data)
  }
  result <- do.call("rbind", result_list)
  names(result) <- c("x", "y", "class")
  return(result)
}