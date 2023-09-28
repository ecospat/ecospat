## Additional functions to the Biomod2 package 13/03/13 Frank Breiner
## Not available anymore until future possible updates
## Rd files
# \name{ecospat.makeDataFrame}
# 
# \alias{ecospat.makeDataFrame}
# 
# \title{Make Data Frame}
# 
# \description{Create a biomod2-compatible dataframe. The function also enables to remove duplicate presences within a pixel and to set a minimum distance between presence points to avoid autocorrelation. Data from GBIF can be added.}
# 
# \usage{ecospat.makeDataFrame (spec.list, expl.var, use.gbif=FALSE, precision=NULL,
#                               year=NULL, remdups=TRUE, mindist=NULL, n=1000, type='random', PApoint=NULL,
#                               ext=expl.var, tryf=5)}
# 
# \arguments{
#   \item{spec.list}{Data.frame or Character. The species occurrence information must be a data.frame in the form:
#       \'x-coordinates\' , \'y-coordinates\' and \'species name\' (in the same projection/coordinate system as expl.var!).}
#   \item{expl.var}{a RasterStack object of the environmental layers.}
#   \item{use.gbif}{Logical. If TRUE presence data from GBIF will be added. It is also possible to use GBIF data only.
#    Default: FALSE. See ?gbif {dismo} for more information. Settings: geo=TRUE, removeZeros=TRUE,
#    all sub-taxa will be used.
#    \'species name\' in spec.list must be in the form: \'genus species\',
#    \'genus_species\' or \'genus.species\'. If there is no species information available on GBIF an error is returned.
#    Try to change species name (maybe there is a synonym) or switch use.gbif off.}
#   \item{precision}{Numeric. Use precision if use.gbif = TRUE to set a minimum precision of
#    the presences which should be added. For precision = 1000
#    e.g. only presences with precision of at least 1000 meter will be added from GBIF.
#    When precision = NULL all presences from GBIF will be used, also presences
#    where precision information is NA.}
#   \item{year}{Numeric. Latest year of the collected gbif occurrences. If year=1960 only occurrences which were collected since 1960 are used.}
#   \item{remdups}{Logical. If TRUE (Default) duplicated presences within a raster pixel
#    will be removed. You will get only one presence per pixel.}
#   \item{mindist}{Numeric. You can set a minimum distance between presence points to avoid
#    autocorrelation. nndist {spatstat} is used to calculate the nearest neighbour (nn)
#    for each point. From the pair of the minimum nn, the point is removed of which the second
#    neighbour is closer. Unit is the same as expl.var.}
#   \item{n}{number of Pseudo-Absences. Default 1000.}
#   \item{type}{sampling dessign for selecting Pseudo-Absences. If \'random\' (default) background points are selected with the 
#    function randomPoints {dismo}. When selecting another sampling type (\'regular\', \'stratified\', \'nonaligned\',
#    \'hexagonal\', \'clustered\' or \'Fibonacci\') spsample {sp} will be used. This can immensely increase computation time and RAM usage if ext is a raster,
#    especially for big raster layers because it must be converted into a  \'SpatialPolygonsDataFrame\' first.}
#   \item{PApoint}{data.frame or SpatialPoints. You can use your own set of Pseudo-Absences
#    instead of generating new PAs. Two columns with \'x\' and \'y\' in the same
#    projection/coordinate system as expl.var!}
#   \item{ext}{a Spatial Object or Raster object. Extent from which PAs should be selected from (Default is expl.var).} 
#   \item{tryf}{numeric > 1. Number of trials to create the requested Pseudo-Absences after removing NA points (if type='random'). See ?randomPoints {dismo}}
# }
# 
# \details{If you use a raster stack as explanatory variable and you want to model many species in a loop with Biomod, formating data will last very long as presences and PA's have to be extracted over and over again from the raster stack. To save computation time, it is better to convert the presences and PAs to a data.frame first.}
# 
# \value{A data.frame object which can be used for modeling with the Biomod package.}
# 
# \author{Frank Breiner \email{frank.breiner@unil.ch}}
# 
# \examples{
#   \donttest{
#     
#     library(dismo)
#     
#     files <- list.files(path=paste(system.file(package="dismo"),
#                                    '/ex', sep=''), pattern='grd', full.names=TRUE )
#     predictors <- raster::stack(files[c(9,1:8)])   #file 9 has more NA values than
#     # the other files, this is why we choose it as the first layer (see ?randomPoints)
#     
#     file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
#     bradypus <- read.table(file, header=TRUE, sep=',')[,c(2,3,1)]
#     head(bradypus)
#     
#     random.spec <- cbind(as.data.frame(randomPoints(predictors,50,extf=1)),species="randomSpec")
#     colnames(random.spec)[1:2] <- c("lon","lat")
#     
#     spec.list <- rbind(bradypus, random.spec)
#     
#     df <- ecospat.makeDataFrame(spec.list, expl.var=predictors, n=5000)
#     head(df)
#     
#     plot(predictors[[1]])
#     points(df[df$Bradypus.variegatus==1, c('x','y')])
#     points(df[df$randomSpec==1, c('x','y')], col="red")
#     
#   }
# }

##############################################

# ecospat.makeDataFrame <- function(spec.list, expl.var, use.gbif = FALSE, precision = NULL,
#                                   year = NULL, remdups = TRUE, mindist = NULL, n = 1000, type = "random", PApoint = NULL,
#                                   ext = expl.var, tryf = 5) {
#   
#   requireNamespace("rgdal")
#   
#   ## Species data
#   if (!is.character(spec.list)) {
#     colnames(spec.list)[1:3] <- c("x", "y", "Spec")
#   }
#   
#   ###sample the Pseudo-Absences
#   
#   if (is.null(PApoint)) {
#     if (type == "random") {
#       n <- n * 1.5
#       PApoint <- data.frame(dismo::randomPoints(mask = expl.var, n, if (!is.character(spec.list)) {
#         p = spec.list[, -3]
#       }, tryf, ext = ext))
#       test <- extract(expl.var, PApoint[, 1:2])
#       n <- (n/1.5)
#       PApoint <- PApoint[complete.cases(test), ][1:n, ]
#     } else {
#       if (inherits(ext,c("RasterStack","RasterLayer"))) {
#         ext <- rasterToPolygons(ext[[1]])
#       }
#       PApoint2 <- data.frame(spsample(ext, n, type))
#     }
#   }
#   PApoint$cell.id <- cellFromXY(expl.var, PApoint)
#   PApoint <- PApoint[!duplicated(PApoint$cell.id), ]
#   PApoint$PA <- 1
#   colnames(PApoint) <- c("x", "y", colnames(PApoint[3:4]))
#   
#   
#   
#   if (!is.character(spec.list)) {
#     # spec.list$Spec <- gsub('\\.', '_', spec.list$Spec) spec.list$Spec <- gsub('
#     # ', '_', spec.list$Spec)
#     spec.list$Spec <- gsub("_", "\\.", spec.list$Spec)
#     spec.list$Spec <- gsub(" ", "\\.", spec.list$Spec)
#     specnames <- levels(factor(spec.list$Spec))
#     
#   } else {
#     specnames <- spec.list
#     specnames <- gsub("_", "\\.", specnames)
#     specnames <- gsub(" ", "\\.", specnames)
#     spec.list <- specnames
#   }
#   
#   
#   ### Beginn Species Loop
#   for (i in 1:length(specnames)) {
#     # plot(i)
#     if (!is.character(spec.list)) {
#       spec.x <- spec.list[, 3] == specnames[i]
#       spec.x <- subset(spec.list, spec.x)
#       spec.data <- data.frame(spec.x[, 3:ncol(spec.x)])
#       colnames(spec.data) <- colnames(spec.x)[3:ncol(spec.x)]
#       spec.x <- spec.x[, c(1:2)]
#     }
#     
#     ### use GBIF data additional to your own
#     if (use.gbif) {
#       sp.name <- gsub(" ", ".", specnames)
#       sp.name <- gsub("_", ".", sp.name)
#       sp.name <- strsplit(sp.name, "\\.")
#       sp.name <- unlist(sp.name)
#       sp.name <- matrix(sp.name, ncol = 2, byrow = TRUE)
#       sp.name <- cbind(specnames, sp.name)
#       
#       spec.gbif <- dismo::gbif(sp.name[i, 2], paste(sp.name[i, 3], "*", sep = ""),
#                                ext = ext, sp = TRUE, geo = TRUE, removeZeros = FALSE, download = TRUE)
#       colnames(spec.gbif@data) <- gsub("species", "Spec", colnames(spec.gbif@data))
#       spec.gbif$coordUncertaintyM <- as.numeric(spec.gbif$coordUncertaintyM)
#       spec.gbif@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#       
#       # set minimum precision of GBIF data
#       if (!is.null(precision)) {
#         prec.na <- is.na(spec.gbif$coordUncertaintyM)
#         if (sum(is.na(spec.gbif$coordUncertaintyM)) == nrow(spec.gbif)) {
#           warning(cat(paste("\n##########\n", "There are no precise presences available at GBIF for species", sp.name[i],", no information about precision.",
#                             "\n##########\n")))
#         }
#         
#         spec.gbif <- spec.gbif[!prec.na, ]
#         minprec <- spec.gbif$coordUncertaintyM <= precision
#         if (sum(minprec) == 0) {
#           warning(cat(paste("\n##########\n", "There are no presences available at GBIF with the selected precision.",
#                             "\n##########\n")))
#         }
#         spec.gbif <- spec.gbif[minprec, ]
#       }
#       # set minimum year of latest Date collected
#       if (!is.null(year)) {
#         spec.gbif <- spec.gbif[!is.na(spec.gbif$latestDateCollected), ]
#         spec.gbif <- spec.gbif[as.numeric(matrix(unlist(strsplit(spec.gbif$latestDateCollected,
#                                                                  "-")), ncol = 3, byrow = TRUE)[, 1]) >= year, ]
#       }
#       
#       # Change projection of GBIF data to projection of your Environmental Raster
#       # Layer.
#       spec.gbif <- spTransform(spec.gbif, CRS(projection(expl.var)))
#       spec.gbif <- spec.gbif[, "Spec"]
#       
#       # Use GBIF data only (if you don't have own occurrences)
#       if (!is.character(spec.list)) {
#         spec.x <- SpatialPointsDataFrame(spec.x[, 1:2], data = spec.data,
#                                          proj4string = CRS(projection(expl.var)))
#         spec.x <- maptools::spRbind(spec.x, spec.gbif)
#         spec.x <- data.frame(spec.x)[, c(2, 3, 1)]
#         spec.x <- spec.x[-which(spec.x$x == 0 & spec.x$y == 0), ]
#       } else {
#         spec.x <- data.frame(spec.gbif)[, c("lon", "lat", "Spec")]
#         spec.x <- spec.x[which(spec.x$lon != 0 & spec.x$lat != 0), ]
#         colnames(spec.x) <- c("x", "y", "Spec")
#       }
#     }
#     # End GBIF
#     
#     ### Check for presence replicates and minimum nearest neighbour
#     spec.x$cell.id <- cellFromXY(expl.var[[1]], spec.x[, c(1, 2)])
#     ### remove Presences with NA-values
#     spec.x <- spec.x[!is.na(spec.x$cell.id), ]
#     
#     # remove duplicated presences within one pixel
#     dup <- duplicated(spec.x$cell.id)
#     if (remdups) {
#       spec.x <- spec.x[!dup, ]
#     } else if (sum(dup) > 0) {
#       warning(cat(paste("\n################################\n", specnames[i], "Warning: There are",
#                         sum(dup), "cells with more than one Presence Points", "\n################################")))
#     }
#     
#     # set a minimum distance between your presences to the nearest neighbour
#     if (!is.null(mindist)) {
#       #spec.x$nn <- nndist(spec.x[, c(1, 2)])
#       #spec.x$nn2 <- nndist(spec.x[, c(1, 2)], k = 2)
#       spec.x$nn <- nabor::knn(spec.x[, c(1, 2)],k=2)$nn.dists[,2]
#       spec.x$nn2 <- nabor::knn(spec.x[, c(1, 2)],k=3)$nn.dists[,3]     
#       
#       spec.x$id <- 1:nrow(spec.x)
#       del.id = 0
#       
#       repeat {
#         if (min(spec.x$nn) < mindist) {
#           mini <- spec.x$nn == min(spec.x$nn)
#           mini2 <- min(spec.x[mini, ]$nn2)
#           del <- spec.x[mini, ]$nn2 == mini2
#           
#           if (sum(del) > 1) {
#             del.id <- spec.x[mini, ][1, ]$id
#           } else del.id <- spec.x[mini, ][del, ]$id
#         }
#         
#         spec.x <- spec.x[!spec.x$id == del.id, ]
#         #spec.x$nn <- nndist(spec.x[, c(1, 2)])
#         #spec.x$nn2 <- nndist(spec.x[, c(1, 2)], k = 2)
#         spec.x$nn <- nabor::knn(spec.x[, c(1, 2)],k=2)$nn.dists[,2]
#         spec.x$nn2 <- nabor::knn(spec.x[, c(1, 2)],k=3)$nn.dists[,3] 
#         
#         if (sum(spec.x$nn < mindist) == 0) {
#           break
#         }
#       }
#     }
#     ### End check for replicates
#     
#     # Print warning massages and summary
#     
#     message(cat(paste("\n\n\n################################\n", specnames[i])))
#     if (use.gbif == TRUE) {
#       message(cat(paste("\n##########\n", "Occurrence data of following species were added from GBIF:","\n")))
#       message(cat(paste(levels(as.factor(spec.x[, 3])), "\n")))
#     }
#     
#     message(cat(paste("\n##########\n", "Dataframe created with ", nrow(spec.x), " Presence Points and ",
#                       sum(PApoint$PA, na.rm = TRUE), " Pseudo-Absence Points.")))
#     
#     
#     if (nrow(spec.x) < 20) {
#       warning(cat(paste("\n##########\n", "Warning: You only have", nrow(spec.x), " Presence Points for modeling.")))
#     }
#     
#     
#     if (nrow(spec.x) < nlayers(expl.var) * 10) {
#       warning(cat(paste("\n##########\n", "Warning: Number of presence points is less than 10 x number of predictors. Be aware of overparametrization. You only have",
#                         nrow(spec.x), "Presences but", nlayers(expl.var), "predictors.")))
#     }
#     message(cat(paste("\n################################\n")))
#     
#     ## Combine Pseudo-Absences and Species data to one data.frame
#     
#     PApoint$specnames <- 0
#     
#     if (i > 1 & sum(is.na(PApoint$PA)) > 0) {
#       PApoint$specnames[(sum(PApoint$PA, na.rm = TRUE) + 1):nrow(PApoint)] <- NA
#     }
#     colnames(PApoint) <- c(colnames(PApoint[1:ncol(PApoint) - 1]), specnames[i])
#     
#     for (j in 1:nrow(spec.x)) {
#       if (sum(PApoint$cell.id == spec.x$cell.id[j]) > 0) {
#         PApoint[which(PApoint$cell.id == spec.x$cell.id[j]), ncol(PApoint)] <- 1
#       } else {
#         add.spec <- data.frame(c(spec.x[j, c("x", "y", "cell.id")], rep(NA,
#                                                                         i), 1))
#         colnames(add.spec) <- colnames(PApoint)
#         PApoint <- rbind(PApoint, add.spec)
#       }
#     }
#     
#   }
#   # End Species Loop. Data Frame created
#   
#   mydataframe <- PApoint
#   
#   expl <- extract(expl.var, mydataframe[, 1:2])
#   mydataframe <- cbind(mydataframe, expl)
#   
#   ### Warning for multicollinearity if |r|>= 0.7
#   if(nlayers(expl.var)>1){
#     if (sum(abs(as.dist(cor(mydataframe[, (ncol(mydataframe) - nlayers(expl.var) +
#                                            1):(ncol(mydataframe))]))) >= 0.7) > 0) {
#       warning(cat(paste("\n\n\n################################\n", "Warning: There are",
#                         sum(abs(as.dist(cor(mydataframe[, (ncol(mydataframe) - nlayers(expl.var) +
#                                                              1):(ncol(mydataframe))]))) >= 0.7), " predictor variable pairs with a correlation coefficients of |r| > 0.7. Be aware of collinearity!",
#                         "\n################################\n")))
#       print(abs(as.dist(round(cor(mydataframe[, (ncol(mydataframe) - nlayers(expl.var) +
#                                                    1):(ncol(mydataframe))]), 2))))
#     }}
#   
#   return(mydataframe)
# }
