## Written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
## University of Lausanne. Switzerland. October 09.
##
## DESCRIPTION
##
## functions to perform measures of niche overlap and niche equivalency/similarity tests as described in Broennimann et al. (submitted)
## 
## list of functions:
##
## grid.clim(glob,glob1,sp,R) 
## use the scores of an ordination (or SDM predictions) and create a grid z of RxR pixels 
## (or a vector of R pixels when using scores of dimension 1 or SDM predictions) with occurrence densities
## Only scores of one, or two dimensions can be used 
## sp= scores for the occurrences of the species in the ordination, glob = scores for the whole studies areas, glob 1 = scores for the range of sp 
##
## niche.overlap(z1,z2,cor)
## calculate the overlap metrics D and I (see Warren et al 2008) based on two species occurrence density grids z1 and z2 created by grid.clim
## cor=T correct occurrence densities of each species by the prevalence of the environments in their range
##
## niche.equivalency.test(z1,z2,rep)
## runs niche equivalency test(see Warren et al 2008) based on two species occurrence density grids
## compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim and z2.sim.
## z1.sim and z2.sim are built from random reallocations of occurences of z1 and z2
## rep is the number of iterations
##
## niche.similarity.test(z1,z2,rep)
## runs niche similarity test(see Warren et al 2008) based on two species occurrence density grids
## compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches (z2.sim) in available in the range of z2 (z2$Z) 
## z2.sim have the same patterns as z2 but their center are randomly translatated in the availabe z2$Z space and weighted by z2$Z densities
## rep is the number of iterations
##
## plot.niche(z,title,name.axis1,name.axis2)
## plot a niche z created by grid.clim. title,name.axis1 and name.axis2 are strings for the legend of the plot
##
## plot.contrib(contrib,eigen)
## plot the contribution of the initial variables to the analysis. Typically the eigen vectors and eigen values in ordinations
##
## plot.overlap.test(x,type,title)
## plot an histogram of observed and randomly simulated overlaps, with p-values of equivalency and similarity tests. 
## x must be an object created by niche.similarity.test or niche.equivalency.test.
## type is either "D" or "I". title is the title of the plot

##################################################################################################

ecospat.niche.overlap <- function(z1, z2, cor) {

  # z1 = species 1 occurrence density grid created by grid.clim
  # z2 = species 2 occurrence density
  # grid created by grid.clim cor= TRUE correct occurrence densities of each species by the prevalence of the environments in their range

  l <- list()

  if (cor == FALSE) {
    p1 <- as.matrix(z1$z.uncor)/sum(as.matrix(z1$z.uncor))  # rescale occurence densities so that the sum of densities is the same for both species
    p2 <- as.matrix(z2$z.uncor)/sum(as.matrix(z2$z.uncor))  # rescale occurence densities so that the sum of densities is the same for both species
  }

  if (cor == TRUE) {
    p1 <- as.matrix(z1$z.cor)/sum(as.matrix(z1$z.cor))  # rescale occurence densities so that the sum of densities is the same for both species
    p2 <- as.matrix(z2$z.cor)/sum(as.matrix(z2$z.cor))  # rescale occurence densities so that the sum of densities is the same for both species
  }

  D <- 1 - (0.5 * (sum(abs(p1 - p2))))  # overlap metric D
  H <- sqrt(sum((sqrt(p1) - sqrt(p2))^2))
  I <- 1 - (H^2)/2  # corrected overlap metric I http://enmtools.blogspot.com.au/2010_09_01_archive.html
  l$D <- D
  l$I <- I
  return(l)
}

################################################################################################## internal function to generate random distribution followint the niche equivalency approach
overlap.eq.gen <- function(repi, z1, z2) {
  if (is.null(z1$y)) {
    # overlap on one axis

    occ.pool <- c(z1$sp, z2$sp)  # pool of random occurrences
    rand.row <- sample(1:length(occ.pool), length(z1$sp))  # random reallocation of occurrences to datasets
    sp1.sim <- occ.pool[rand.row]
    sp2.sim <- occ.pool[-rand.row]
  }

  if (!is.null(z1$y)) {
    # overlap on two axes

    occ.pool <- rbind(z1$sp, z2$sp)  # pool of random occurrences
    row.names(occ.pool)<-c()  # remove the row names
    rand.row <- sample(1:nrow(occ.pool), nrow(z1$sp))  # random reallocation of occurrences to datasets
    sp1.sim <- occ.pool[rand.row, ]
    sp2.sim <- occ.pool[-rand.row, ]
  }

  z1.sim <- ecospat.grid.clim.dyn(z1$glob, z1$glob1, data.frame(sp1.sim), R = length(z1$x))  # gridding
  z2.sim <- ecospat.grid.clim.dyn(z2$glob, z2$glob1, data.frame(sp2.sim), R = length(z2$x))

  o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = TRUE)  # overlap between random and observed niches
  sim.o.D <- o.i$D  # storage of overlaps
  sim.o.I <- o.i$I
  return(c(sim.o.D, sim.o.I))
}


ecospat.niche.equivalency.test <- function(z1, z2, rep, alternative = "greater", ncores=1) {

  R <- length(z1$x)
  l <- list()

  obs.o <- ecospat.niche.overlap(z1, z2, cor = TRUE)  #observed niche overlap

  if (ncores == 1){
    sim.o <- as.data.frame(matrix(unlist(lapply(1:rep, overlap.eq.gen, z1, z2)), byrow = TRUE,
      ncol = 2))  #simulate random overlap
    }else{
  #number of cores attributed for the permutation test
  cl <- makeCluster(ncores)  #open a cluster for parallelization
  invisible(clusterEvalQ(cl, library(ecospat)))  #import the internal function into the cluster
  sim.o <- as.data.frame(matrix(unlist(parLapply(cl, 1:rep, overlap.eq.gen, z1, z2)), byrow = TRUE,
    ncol = 2))  #simulate random overlap
  stopCluster(cl)  #shutdown the cluster
    }
  colnames(sim.o) <- c("D", "I")
  l$sim <- sim.o  # storage
  l$obs <- obs.o  # storage

  if (alternative == "greater") {
    l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values alternative hypothesis = greater -> test for niche conservatism/convergence
    l$p.I <- (sum(sim.o$I >= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of p-values alternative hypothesis = greater -> test for niche conservatism/convergence
  }
  if (alternative == "lower") {
    l$p.D <- (sum(sim.o$D <= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values alternative hypothesis = lower -> test for niche divergence
    l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of p-values alternative hypothesis = lower -> test for niche divergence
  }

  return(l)
}

##################################################################################################

#### internal function to generate random distribution following the niche similarity approach
overlap.sim.gen <- function(repi, z1, z2, rand.type = rand.type) {
  R1 <- length(z1$x)
  R2 <- length(z2$x)
  if (is.null(z1$y) & is.null(z2$y)) {
    if (rand.type == 1) {
      # if rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly
      # shifted
      center.z1 <- which(z1$z.uncor == 1)  # define the centroid of the observed niche
      Z1 <- z1$Z/max(z1$Z)
      rand.center.z1 <- sample(1:R1, size = 1, replace = FALSE, prob = Z1)  # randomly (weighted by environment prevalence) define the new centroid for the niche
      xshift.z1 <- rand.center.z1 - center.z1  # shift on x axis
      z1.sim <- z1
      z1.sim$z <- rep(0, R1)  # set intial densities to 0
      for (i in 1:length(z1$x)) {
        i.trans.z1 <- i + xshift.z1
        if (i.trans.z1 > R1 | i.trans.z1 < 0)
          (next)()  # densities falling out of the env space are not considered
        z1.sim$z[i.trans.z1] <- z1$z[i]  # shift of pixels
      }
      z1.sim$z <- (z1$Z != 0) * 1 * z1.sim$z  # remove densities out of existing environments
      z1.sim$z.cor <- (z1.sim$z/z1$Z)/max((z1.sim$z/z1$Z), na.rm = TRUE)  #transform densities into occupancies
      z1.sim$z.cor[which(is.na(z1.sim$z.cor))] <- 0
      z1.sim$z.uncor <- z1.sim$z/max(z1.sim$z, na.rm = TRUE)
      z1.sim$z.uncor[which(is.na(z1.sim$z.uncor))] <- 0
    }

    center.z2 <- which(z2$z.uncor == 1)  # define the centroid of the observed niche
    Z2 <- z2$Z/max(z2$Z)
    rand.center.z2 <- sample(1:R2, size = 1, replace = FALSE, prob = Z2)  # randomly (weighted by environment prevalence) define the new centroid for the niche

    xshift.z2 <- rand.center.z2 - center.z2  # shift on x axis
    z2.sim <- z2
    z2.sim$z <- rep(0, R2)  # set intial densities to 0
    for (i in 1:length(z2$x)) {
      i.trans.z2 <- i + xshift.z2
      if (i.trans.z2 > R2 | i.trans.z2 < 0)
        (next)()  # densities falling out of the env space are not considered
      z2.sim$z[i.trans.z2] <- z2$z[i]  # shift of pixels
    }
    z2.sim$z <- (z2$Z != 0) * 1 * z2.sim$z  # remove densities out of existing environments
    z2.sim$z.cor <- (z2.sim$z/z2$Z)/max((z2.sim$z/z2$Z), na.rm = TRUE)  #transform densities into occupancies
    z2.sim$z.cor[which(is.na(z2.sim$z.cor))] <- 0
    z2.sim$z.uncor <- z2.sim$z/max(z2.sim$z, na.rm = TRUE)
    z2.sim$z.uncor[which(is.na(z2.sim$z.uncor))] <- 0
  }

  if (!is.null(z2$y) & !is.null(z1$y)) {
    if (rand.type == 1) {
      # if rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly
      # shifted
      centroid.z1 <- which(z1$z.uncor == 1, arr.ind = TRUE)[1, ]  # define the centroid of the observed niche
      Z1 <- z1$Z/max(z1$Z)
      rand.centroids.z1 <- which(Z1 > 0, arr.ind = TRUE)  # all pixels with existing environments in the study area
      weight.z1 <- Z1[Z1 > 0]
      rand.centroid.z1 <- rand.centroids.z1[sample(1:nrow(rand.centroids.z1), size = 1, replace = FALSE,
        prob = weight.z1), ]  # randomly (weighted by environment prevalence) define the new centroid for the niche
      xshift.z1 <- rand.centroid.z1[1] - centroid.z1[1]  # shift on x axis
      yshift.z1 <- rand.centroid.z1[2] - centroid.z1[2]  # shift on y axis
      z1.sim <- z1
      z1.sim$z <- matrix(rep(0, R1 * R1), ncol = R1, nrow = R1)  # set intial densities to 0
      for (i in 1:R1) {
        for (j in 1:R1) {
          i.trans.z1 <- i + xshift.z1
          j.trans.z1 <- j + yshift.z1
          if (i.trans.z1 > R1 | i.trans.z1 < 0)
          (next)()  # densities falling out of the env space are not considered
          if (j.trans.z1 > R1 | j.trans.z1 < 0)
          (next)()
          z1.sim$z[i.trans.z1, j.trans.z1] <- z1$z[i, j]  # shift of pixels
        }
      }
      z1.sim$z <- (z1$Z != 0) * 1 * z1.sim$z  # remove densities out of existing environments
      z1.sim$z.cor <- (z1.sim$z/z1$Z)/max((z1.sim$z/z1$Z), na.rm = TRUE)  #transform densities into occupancies
      z1.sim$z.cor[which(is.na(z1.sim$z.cor))] <- 0
      z1.sim$z.uncor <- z1.sim$z/max(z1.sim$z, na.rm = TRUE)
      z1.sim$z.uncor[which(is.na(z1.sim$z.uncor))] <- 0
    }
    centroid.z2 <- which(z2$z.uncor == 1, arr.ind = TRUE)[1, ]  # define the centroid of the observed niche
    Z2 <- z2$Z/max(z2$Z)
    rand.centroids.z2 <- which(Z2 > 0, arr.ind = TRUE)  # all pixels with existing environments in the study area
    weight.z2 <- Z2[Z2 > 0]
    rand.centroid.z2 <- rand.centroids.z2[sample(1:nrow(rand.centroids.z2), size = 1, replace = FALSE,
      prob = weight.z2), ]  # randomly (weighted by environment prevalence) define the new centroid for the niche
    xshift.z2 <- rand.centroid.z2[1] - centroid.z2[1]  # shift on x axis
    yshift.z2 <- rand.centroid.z2[2] - centroid.z2[2]  # shift on y axis
    z2.sim <- z2
    z2.sim$z <- matrix(rep(0, R2 * R2), ncol = R2, nrow = R2)  # set intial densities to 0
    for (i in 1:R2) {
      for (j in 1:R2) {
        i.trans.z2 <- i + xshift.z2
        j.trans.z2 <- j + yshift.z2
        if (i.trans.z2 > R2 | i.trans.z2 < 0)
          (next)()  # densities falling out of the env space are not considered
        if (j.trans.z2 > R2 | j.trans.z2 < 0)
          (next)()
        z2.sim$z[i.trans.z2, j.trans.z2] <- z2$z[i, j]  # shift of pixels
      }
    }
    z2.sim$z <- (z2$Z != 0) * 1 * z2.sim$z  # remove densities out of existing environments
    z2.sim$z.cor <- (z2.sim$z/z2$Z)/max((z2.sim$z/z2$Z), na.rm = TRUE)  #transform densities into occupancies
    z2.sim$z.cor[which(is.na(z2.sim$z.cor))] <- 0
    z2.sim$z.uncor <- z2.sim$z/max(z2.sim$z, na.rm = TRUE)
    z2.sim$z.uncor[which(is.na(z2.sim$z.uncor))] <- 0
  }

  if (rand.type == 1) {
    o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = TRUE)
  }
  if (rand.type == 2)
    {
      o.i <- ecospat.niche.overlap(z1, z2.sim, cor = TRUE)
    }  # overlap between random and observed niches
  sim.o.D <- o.i$D  # storage of overlaps
  sim.o.I <- o.i$I
  return(c(sim.o.D, sim.o.I))
}

ecospat.niche.similarity.test <- function(z1, z2, rep, alternative = "greater", rand.type = 1, ncores = 1) {

  R <- length(z1$x)
  l <- list()
  obs.o <- ecospat.niche.overlap(z1, z2, cor = TRUE)  #observed niche overlap
  z1$z.uncor <- as.matrix(z1$z.uncor)
  z1$Z <- as.matrix(z1$Z)
  z1$z <- as.matrix(z1$z)
  z2$z.uncor <- as.matrix(z2$z.uncor)
  z2$Z <- as.matrix(z2$Z)
  z2$z <- as.matrix(z2$z)

  if (ncores==1) {
    sim.o <- as.data.frame(matrix(unlist(lapply(1:rep, overlap.sim.gen, z1, z2, rand.type = rand.type)),
                                  byrow = TRUE, ncol = 2))  #simulate random overlap  
    } else {
  cl <- makeCluster(ncores)  #open a cluster for parallelization
  invisible(clusterEvalQ(cl, library(ecospat)))  #import the internal function into the cluster
  sim.o <- as.data.frame(matrix(unlist(parLapply(cl, 1:rep, overlap.sim.gen, z1, z2, rand.type = rand.type)),
    byrow = TRUE, ncol = 2))  #simulate random overlap
  stopCluster(cl)  #shutdown the cluster
    }
  colnames(sim.o) <- c("D", "I")
  l$sim <- sim.o  # storage
  l$obs <- obs.o  # storage

  if (alternative == "greater") {
    l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values alternative hypothesis = greater -> test for niche conservatism/convergence
    l$p.I <- (sum(sim.o$I >= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of p-values alternative hypothesis = greater -> test for niche conservatism/convergence
  }
  if (alternative == "lower") {
    l$p.D <- (sum(sim.o$D <= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values alternative hypothesis = lower -> test for niche divergence
    l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of p-values alternative hypothesis = lower -> test for niche divergence
  }

  return(l)
}

##################################################################################################

ecospat.plot.niche <- function(z, title = "", name.axis1 = "Axis 1", name.axis2 = "Axis 2", cor = FALSE) {
  if (is.null(z$y)) {
    R <- length(z$x)
    x <- z$x
    xx <- sort(rep(1:length(x), 2))
    if (cor == FALSE)
      y1 <- z$z.uncor/max(z$z.uncor)
    if (cor == TRUE)
      y1 <- z$z.cor/max(z$z.cor)
    Y1 <- z$Z/max(z$Z)
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 2)]
    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = "grey")
    lines(x[xx], c(0, Y1[YY1], 0, 0))
  }
  if (!is.null(z$y)) {
    if (cor == FALSE)
      image(x=z$x,y=z$y,z=t(as.matrix(z$z.uncor))[,nrow(as.matrix(z$z.uncor)):1], col = gray(100:0/100), zlim = c(1e-06, cellStats(z$z.uncor,"max")),
        xlab = name.axis1, ylab = name.axis2)
    if (cor == TRUE)
      image(x=z$x,y=z$y,z=t(as.matrix(z$z.uncor))[,nrow(as.matrix(z$z.uncor)):1], col = gray(100:0/100), zlim = c(1e-06, cellStats(z$z.cor,"max")), 
      xlab = name.axis1, ylab = name.axis2)
    z$Z<-t(as.matrix(z$Z))[,nrow(as.matrix(z$Z)):1]
    contour(x=z$x,y=z$y,z$Z, add = TRUE, levels = quantile(z$Z[z$Z > 0], c(0, 0.5)), drawlabels = FALSE,
      lty = c(1, 2))
  }
  title(title)
}



ecospat.plot.contrib <- function(contrib, eigen) {

  if (ncol(contrib) == 1) {
    h <- c(unlist(contrib))
    n <- row.names(contrib)
    barplot(h, space = 0, names.arg = n)
    title(main = "variable contribution")
  }
  if (ncol(contrib) == 2) {
    s.corcircle(contrib[, 1:2]/max(abs(contrib[, 1:2])), grid = FALSE)
    title(main = "correlation circle", sub = paste("axis1 = ", round(eigen[1]/sum(eigen) *
      100, 2), "%", "axis2 = ", round(eigen[2]/sum(eigen) * 100, 2), "%"))
  }
}


ecospat.plot.overlap.test <- function(x, type, title) {
  if (type == "D") {
    obs <- x$obs$D
    sim <- x$sim$D
    p <- x$p.D
  }
  if (type == "I") {
    obs <- x$obs$I
    sim <- x$sim$I
    p <- x$p.I
  }
  r0 <- c(sim, obs)
  l0 <- max(sim) - min(sim)
  w0 <- l0/(log(length(sim), base = 2) + 1)
  xlim0 <- range(r0) + c(-w0, w0)
  h0 <- hist(sim, plot = FALSE, nclass = 10)
  y0 <- max(h0$counts)
  hist(sim, plot = TRUE, nclass = 10, xlim = xlim0, col = grey(0.8), main = title, xlab = type,
    sub = paste("p.value = ", round(p, 5)))
  lines(c(obs, obs), c(y0/2, 0), col = "red")
  points(obs, y0/2, pch = 18, cex = 2, col = "red")
  invisible()
}

##################################################################################################