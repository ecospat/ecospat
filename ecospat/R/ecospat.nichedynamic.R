## Written by Olivier Broennimann and Blaise Petitpierre. Departement of Ecology and Evolution (DEE). 
## University of Lausanne. Switzerland. April 2012.
##
## DESCRIPTION
##
## functions to perform measures of niche overlap and niche equivalency/similarity tests as described in Broennimann et al. (submitted)
## 
## list of functions:
##
## ecospat.kd(x,ext,R = 100,th = 0 ,env.mask = c(),method = 'adehabitat')
## wrapper function to draw smoothed densities in the environmental space using different methods
## x = environmental  scores (along one or two dimension)
## ext = extent along the environmental score. Vector including corresponding to c(xmin, xmax, ymin, ymax). The y component is optional.
## th = percentile threshold used to delineate the niche. 0 means that all occurences are included within the niche. 0.05 means that 95 % of the distribution is included within the niche.
## env.mask = raster of vector containing the environmental background densities
## method = method used to draw the kernel density distribution. By default, 'adehabitat' uses kernelUD from the library adehabitatHR but you can set it to 'ks' to use the the algorithm from the ks library
## 
## grid.clim.dyn(glob,glob1,sp,R,th.sp,th.env) 
## use the scores of an ordination (or SDM predictions) and create a grid z of RxR pixels 
## (or a vector of R pixels when using scores of dimension 1 or SDM predictions) with occurrence densities
## Only scores of one, or two dimensions can be used 
## sp= scores for the occurrences of the species in the ordination, glob = scores for the whole studies areas, glob 1 = scores for the range of sp 
## R= resolution of the grid, th.sp=quantile of species densitie at species occurences used as a threshold to exclude low species density values, 
## th.env=The quantile used to delimit a threshold to exclude low environmental density values of the study area.
## geomask= a geographical mask to delimit the background extent if the analysis takes place in the geographical space. It can be a SpatialPolygon or a raster object. Note that the CRS should be the same as the one used for the points.
## kernel.method = method used to estimate the the kernel density. The initial and original method is 'adehabitat', while 'ks' has been recently implemented for future developments in multidimensional space
## 
## dynamic.index(z1,z2,intersection=NA)
## calculate niche expansion, stability and unfilling
## z1 : gridclim object for the native distribution
## z2 : gridclim object for the invaded range
## intersection : quantile of the environmental density used to remove marginal climates. 
## If intersection = NA, analysis is performed on the whole environmental extent (native and invaded)
## If intersection = 0, analysis is performed at the intersection between native and invaded range
## If intersection = 0.05, analysis is performed at the intersection of the 5th quantile of both native and invaded environmental densities 
## etc...
##
## plot.niche.dyn(z1,z2,quant,title,interest,colz1,colz2,colinter,colZ1=,colZ2=)
## plot niche categories and species density
## z1 : gridclim object for the native distribution
## z2 : gridclim object for the invaded range
## quant : quantile of the environmental density used to delimit marginal climates.
## title : title of the figure
## interest : choose which density to plot. If interest=1 plot native density, if interest=2 plot invasive density
## colz1 : color used to depict unfilling area
## colz2 : color used to depict expansion area
## colinter : color used to depict overlap area
## colZ1 : color used to delimit the native extent
## colZ2 : color used to delimit the invaded extent
##
## ecospat.shift.centroids(sp1,sp2,clim1,clim2)
## draw arrows linking the centroid of the native and inasive distribution (continuous line) and between native and invaded extent (dashed line)
## sp1 : scores of the species native distribution along the the 2 first axes of the PCA
## sp2 : scores of the species invasive distribution along the the 2 first axes of the PCA
## clim1 : scores of the entire native extent along the the 2 first axes of the PCA
## clim2 : scores of the entire invaded extent along the the 2 first axes of the PCA


##################################################################################################

ecospat.kd<-function(x,ext,R = 100,th = 0,env.mask = c(),
                     method = 'adehabitat'){
  if (method == 'adehabitat'){
    if (ncol(x) == 2){
      xr <- data.frame(cbind((x[, 1] - ext[1])/abs(ext[2] - ext[1]), 
                             (x[, 2] -ext[3])/abs(ext[4] - ext[3])))  # data preparation
      mask <- ascgen(SpatialPoints(cbind((0:(R))/R, (0:(R)/R))), 
                     nrcol = R-2, count = FALSE) # data preparation
      x.dens <- kernelUD(SpatialPoints(xr[, 1:2]), h = "href", grid = mask,
                         kern = "bivnorm")  # calculate the density of occurrences in a grid of RxR pixels along the score gradients
      x.dens <- raster(xmn = ext[1], xmx = ext[2], ymn = ext[3], 
                       ymx = ext[4],matrix(x.dens$ud, nrow = R))
      if (!is.null(th)){
        th.value<-quantile(extract(x.dens,x),th)
        x.dens[x.dens<th.value]<-0
      }
      if (!is.null(env.mask)){
        x.dens<-x.dens * env.mask    
      }
    }  else if (ncol(x) == 1){  
      xr <- seq(from = min(ext), to = max(ext), length.out = R)  # breaks on score gradient 1
      x.dens <- density(x[, 1], kernel = "gaussian", from = min(xr), to = max(xr),
                        n = R, cut = 0)  # calculate the density of occurrences in a vector of R pixels along the score gradient
      # using a gaussian kernel density function, with R bins.
      if (!is.null(env.mask)){
        x.dens$y<-x.dens$y * env.mask    
      }
      if (!is.null(th)){
        xr <- sapply(x, findInterval, x.dens$x)
        th.value <- quantile(x.dens$y[xr], th)
        sprm <- which(x.dens$y < th.value)
        x.dens$y[sprm] <- 0  # remove infinitesimally small number generated by kernel density function
      }      
    } 
  }
  
  if (method=='ks'){
    if (ncol(x) == 2){    
      x.dens<-kde(x,xmin=ext[c(1,3)],
                  xmax = ext[c(2,4)],gridsize =c(R,R))
      x.dens <- flip(t(raster(x.dens$estimate)),direction = 'y')
      extent(x.dens)<-c(xmn = ext[1], xmx = ext[2], ymn = ext[3],
                        ymx = ext[4])
      if (!is.null(th)){
        th.value<-quantile(extract(x.dens,x),th)
        x.dens[x.dens<th.value]<-0
      }
      if (!is.null(env.mask)){
        x.dens<-x.dens * env.mask    
      }
    }else if (ncol(x) == 1){
      x.dens<-kde(x,xmin=min(ext),
                  xmax = max(ext),gridsize =c(R,R))
      x.dens$y<-x.dens$estimate
      x.dens$x<-x.dens$eval.points
      if (!is.null(env.mask)){
        x.dens$y<-x.dens$y * env.mask    
      }
      if (!is.null(th)){
        xr <- sapply(x, findInterval, x.dens$x)
        th.value <- quantile(x.dens$y[xr], th)
        sprm <- which(x.dens$y < th.value)
        x.dens$y[sprm] <- 0  # remove infinitesimally small number generated by kernel density function
      }      
    }
    
  }
  
  return (x.dens)
}  


##################################################################################################

ecospat.grid.clim.dyn <- function(glob, glob1, sp, R = 100, th.sp = 0, 
                                  th.env = 0, geomask = NULL, kernel.method = 'adehabitat') {
  if (is.null(kernel.method)|(kernel.method!='ks'& kernel.method!='adehabitat'))
    stop("supply a kernel method ('adehabitat' or 'ks')")
  
  glob <- as.matrix(glob)
  glob1 <- as.matrix(glob1)
  sp <- as.matrix(sp)
  l <- list()
  
  if (ncol(glob) > 2)
    stop("cannot calculate overlap with more than two axes")
  
  if (ncol(glob) == 1) {
    # if scores in one dimension (e.g. LDA,SDM predictions,...)
    xmax <- max(glob[, 1])
    xmin <- min(glob[, 1])
    if (kernel.method == 'ks'){
      glob1.dens<-ecospat.kd(x = glob1,ext = c(xmin,xmax),method = 'ks',th=0)
      sp.dens<-ecospat.kd(x = sp,ext = c(xmin,xmax),method = 'ks',th=0,
                          env.mask = glob1.dens$y>0)
    }else if (kernel.method == 'adehabitat'){
      glob1.dens<-ecospat.kd(x = glob1,ext = c(xmin,xmax),method = 'adehabitat',th=0)
      sp.dens<-ecospat.kd(x = sp,ext = c(xmin,xmax),method = 'adehabitat',th=0,
                          env.mask = glob1.dens$y>0)
    }
    x<-sp.dens$x
    y<-sp.dens$y
    z <- sp.dens$y * nrow(sp)/sum(sp.dens$y)  # rescale density to the number of occurrences in sp, ie. number of occurrence/pixel
    Z <- glob1.dens$y * nrow(glob)/sum(glob1.dens$y)  # rescale density to the number of sites in glob1
    
    z.uncor <- z/max(z)  # rescale between [0:1] for comparison with other species
    z.cor <- z/Z  # correct for environment prevalence
    z.cor[is.na(z.cor)] <- 0  # remove n/0 situations
    z.cor[z.cor == "Inf"] <- 0  # remove 0/0 situations
    z.cor <- z.cor/max(z.cor)  # rescale between [0:1] for comparison with other species
  }
  
  if (ncol(glob) == 2) {
    # if scores in two dimensions (e.g. PCA)
    xmin<-apply(glob,2,min,na.rm=T)
    xmax<-apply(glob,2,max,na.rm=T)
    ext = c(xmin[1],xmax[1],xmin[2],xmax[2])
    
    if (kernel.method == 'ks'){
      glob1.dens<-ecospat.kd(x = glob1,ext = ext,method = 'ks',th=0)
      if (!is.null(geomask)) {
        proj4string(geomask) <- NA
        glob1.dens <- mask(glob1.dens, geomask, updatevalue = 0)  # Geographical mask in the case if the analysis takes place in the geographical space
      }
      sp.dens<-ecospat.kd(x = sp,ext = ext ,method = 'ks',th=0,
                          env.mask = glob1.dens>0)
    }else if (kernel.method == 'adehabitat'){
      glob1.dens<-ecospat.kd(x = glob1,ext = ext,method = 'adehabitat',th=0)
      if (!is.null(geomask)) {
        proj4string(geomask) <- NA
        glob1.dens <- mask(glob1.dens, geomask, updatevalue = 0)  # Geographical mask in the case if the analysis takes place in the geographical space
      }
      sp.dens<-ecospat.kd(x = sp,ext = ext,method = 'adehabitat',th=0,
                          env.mask = glob1.dens>0)
    }
    
    x<-seq(from = ext[1],to = ext[2],length.out = 100)
    y<-seq(from = ext[3],to = ext[4],length.out = 100)
    Z <- glob1.dens * nrow(glob1)/cellStats(glob1.dens, "sum") # rescale density to the number of occurrences in sp, ie. number of occurrence/pixel
    z <- sp.dens * nrow(sp)/cellStats(sp.dens, "sum") # rescale density to the number of occurrences in sp, ie. number of occurrence/pixel
    z.uncor <- z/cellStats(z, "max")
    z.cor <- z/Z  # correct for environment prevalence
    z.cor[is.na(z.cor)] <- 0  # remove n/0 situations
    z.cor <- z.cor/cellStats(z.cor, "max")
  }
  
  w <- z.uncor  # niche envelope
  w[w > 0] <- 1
  l$x <- x
  l$y <- y
  l$z <- z
  l$z.uncor <- z.uncor
  l$z.cor <- z.cor
  l$Z <- Z
  l$glob <- glob
  l$glob1 <- glob1
  l$sp <- sp
  l$w <- w
  
  return(l)
}


##################################################################################################

ecospat.plot.niche.dyn <- function(z1, z2, quant, title = "", name.axis1 = "Axis 1",
  name.axis2 = "Axis 2", interest = 1, colz1 = "#00FF0050", colz2 = "#FF000050",
  colinter = "#0000FF50", colZ1 = "green3", colZ2 = "red3") {

  if (is.null(z1$y)) {
    R <- length(z1$x)
    x <- z1$x
    xx <- sort(rep(1:length(x), 2))

    y1 <- z1$z.uncor/max(z1$z.uncor)
    Y1 <- z1$Z/max(z1$Z)
    if (quant > 0) {
      Y1.quant <- quantile(z1$Z[which(z1$Z > 0)], probs = seq(0, 1, quant))[2]/max(z1$Z)
    } else {
      Y1.quant <- 0
    }
    Y1.quant <- Y1 - Y1.quant
    Y1.quant[Y1.quant < 0] <- 0
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 2)]

    y2 <- z2$z.uncor/max(z2$z.uncor)
    Y2 <- z2$Z/max(z2$Z)
    if (quant > 0) {
      Y2.quant <- quantile(z2$Z[which(z2$Z > 0)], probs = seq(0, 1, quant))[2]/max(z2$Z)
    } else {
      Y2.quant = 0
    }
    Y2.quant <- Y2 - Y2.quant
    Y2.quant[Y2.quant < 0] <- 0
    yy2 <- sort(rep(1:length(y2), 2))[-c(1:2, length(y2) * 2)]
    YY2 <- sort(rep(1:length(Y2), 2))[-c(1:2, length(Y2) * 2)]

    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = colz1, border = 0)
    polygon(x[xx], c(0, y2[yy2], 0, 0), col = colz2, border = 0)
    polygon(x[xx], c(0, apply(cbind(y2[yy2], y1[yy1]), 1, min, na.exclude = TRUE),
      0, 0), col = colinter, border = 0)
    lines(x[xx], c(0, Y2.quant[YY2], 0, 0), col = colZ2, lty = "dashed")
    lines(x[xx], c(0, Y1.quant[YY1], 0, 0), col = colZ1, lty = "dashed")
    lines(x[xx], c(0, Y2[YY2], 0, 0), col = colZ2)
    lines(x[xx], c(0, Y1[YY1], 0, 0), col = colZ1)
    segments(x0 = 0, y0 = 0, x1 = max(x[xx]), y1 = 0, col = "white")
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, col = "white")

    seg.cat <- function(inter, cat, col.unf, col.exp, col.stab) {
      if (inter[3] == 0) {
        my.col = 0
      }
      if (inter[3] == 1) {
        my.col = col.unf
      }
      if (inter[3] == 2) {
        my.col = col.stab
      }
      if (inter[3] == -1) {
        my.col = col.exp
      }
      segments(x0 = inter[1], y0 = -0.01, y1 = -0.01, x1 = inter[2],
        col = my.col, lwd = 4, lty = 2)
    }
    cat <- ecospat.niche.dyn.index(z1, z2, intersection = quant)$dyn
    inter <- cbind(z1$x[-length(z1$x)], z1$x[-1], cat[-1])
    apply(inter, 1, seg.cat, col.unf = "#00FF0050", col.exp = "#FF000050",
      col.stab = "#0000FF50")

  }

  if (!is.null(z1$y)) {
    z <- t(as.matrix(z1$w + 2 * z2$w))[,nrow(as.matrix(z1$z.uncor)):1]
    z1$Z<-t(as.matrix(z1$Z))[,nrow(as.matrix(z1$Z)):1]
    z2$Z<-t(as.matrix(z2$Z))[,nrow(as.matrix(z2$Z)):1]
    if (interest == 1) {
      image(x=z1$x,y=z1$y,z=t(as.matrix(z1$z.uncor))[,nrow(as.matrix(z1$z.uncor)):1], col = gray(100:0/100), zlim = c(1e-05, cellStats(z1$z.uncor,"max")), xlab = name.axis1, ylab = name.axis2)
      image(x=z1$x,y=z1$y,z=z, col = c("#FFFFFF00", colz1, colz2, colinter), add = TRUE)
    }
    if (interest == 2) {
      image(x=z2$x,y=z2$y,z=t(as.matrix(z2$z.uncor))[,nrow(as.matrix(z2$z.uncor)):1], col = gray(100:0/100), zlim = c(1e-05, cellStats(z2$z.uncor,"max")), xlab = name.axis1, ylab = name.axis2)
      image(x=z2$x,y=z2$y,z=z, col = c("#FFFFFF00", colz1, colz2, colinter), add = TRUE)
    }
    title(title)
    contour(x=z1$x,y=z1$y,z1$Z, add = TRUE, levels = quantile(z1$Z[z1$Z > 0], c(0, quant)),
      drawlabels = FALSE, lty = c(1, 2), col = colZ1)
    contour(x=z2$x,y=z2$y,z2$Z, add = TRUE, levels = quantile(z2$Z[z2$Z > 0], c(0, quant)),
      drawlabels = FALSE, lty = c(1, 2), col = colZ2)
  }
}

##################################################################################################

ecospat.shift.centroids <- function(sp1, sp2, clim1, clim2, col = "red") {

  if (ncol(as.matrix(sp1)) == 2) {
    arrows(median(sp1[, 1]), median(sp1[, 2]), median(sp2[, 1]), median(sp2[,
      2]), col = "red", lwd = 2, length = 0.1)
    arrows(median(clim1[, 1]), median(clim1[, 2]), median(clim2[, 1]),
      median(clim2[, 2]), lty = "11", col = col, lwd = 2, length = 0.1)
  } else {
    arrows(median(sp1), 0.025, median(sp2), 0.025, col = "red", lwd = 2,
      length = 0.1)
    arrows(median(clim1), -0.025, median(clim2), -0.025, lty = "11", col = col,
      lwd = 2, length = 0.1)
  }
}

##################################################################################################

ecospat.niche.dyn.index <- function(z1, z2, intersection = NA) {
  rotate <- function(x) t(apply(x, 2, rev))
  w1 <- as.matrix(z1$w)  # native environmental distribution mask
  w2 <- as.matrix(z2$w)  # invaded environmental distribution mask
  glob1 <- as.matrix(z1$Z)  # Native environmental extent densities
  glob2 <- as.matrix(z2$Z)  # Invaded environmental extent densities
  if (!is.na(intersection)) {
    if (intersection == 0) {
      glob1[glob1 > 0] <- 1  # Native environmental extent mask
      glob2[glob2 > 0] <- 1  # Invaded environmental extent mask
    } else {
      quant.val <- quantile(glob1[glob1 > 0], probs = seq(0, 1, intersection))[2]  # threshold do delimit native environmental mask
      glob1[glob1[] <= quant.val] <- 0
      glob1[glob1[] > quant.val] <- 1  #  native environmental mask
      quant.val <- quantile(glob2[glob2 > 0], probs = seq(0, 1, intersection))[2]  # threshold do delimit invaded environmental mask
      glob2[glob2[] <= quant.val] <- 0
      glob2[glob2[] > quant.val] <- 1  #  invaded environmental mask
    }
    glob <- glob1 * glob2  # delimitation of the intersection between the native and invaded extents
    w1 <- w1 * glob  # Environmental native distribution at the intersection
    w2 <- w2 * glob  # Environmental invasive distribution at the intersection
  }
  z.exp.cat <- (w1 + 2 * w2)/2
  z.exp.cat[z.exp.cat != 1] <- 0  #categorizing expansion pixels
  z.stable.cat <- (w1 + 2 * w2)/3
  z.stable.cat[z.stable.cat != 1] <- 0  #categorizing stable pixels
  z.res.cat <- w1 + 2 * w2
  z.res.cat[z.res.cat != 1] <- 0  #categorizing restriction pixels
  obs.exp <- as.matrix(z2$z.uncor) * as.matrix(z.exp.cat)  #density correction
  obs.stab <- as.matrix(z2$z.uncor) * as.matrix(z.stable.cat)  #density correction
  obs.res <- as.matrix(z1$z.uncor) * as.matrix(z.res.cat)  #density correction

  dyn <- (-1 * z.exp.cat) + (2 * z.stable.cat) + z.res.cat
  if (ncol(w1) == 2)
    {
      dyn <- raster(dyn)
    }  # draw matrix with 3 categories of niche dynamic
  expansion.index.w <- sum(obs.exp)/sum(obs.stab + obs.exp)  # expansion
  stability.index.w <- sum(obs.stab)/sum(obs.stab + obs.exp)  # stability
  restriction.index.w <- sum(obs.res)/sum(obs.res + (z.stable.cat * as.matrix(z1$z.uncor)))  #unfilling
  part <- list()
  part$dyn <- rotate(dyn)
  part$dynamic.index.w <- c(expansion.index.w, stability.index.w, restriction.index.w)
  names(part$dynamic.index.w) <- c("expansion", "stability", "unfilling")
  return(part)
}
