\name{ecospat.margin}

\alias{ecospat.margin}

\title{Delineation of the distribution's margin and its uncertainty}

\description{Estimate the margin of the distribution in the bi-dimensional environmental space based on bootstrapped kernel density estimation the percentage of the distribution included within the margin.}

\usage{ecospat.margin (z, th.quant, kern.method, bootstrap, bootstrap.rep, bootstrap.ncore)}

\arguments{
  \item{z}{object (list) resulting from the function \code{ecospat.grid.clim.dyn}}
  \item{th.quant}{The quantile (between 0 and 1) used to delimit a threshold to exclude low species density values (see details)}
  \item{kern.method}{Method used to draw the kernel density estimation. it can be "adehabitat" (default) or "ks"}
  \item{bootstrap}{Boolean. If \code{TRUE}, a confidence interval based on bootstrap is estimated for the margin of the distribution}
  \item{bootstrap.rep}{Number of resamplings if the boostrap is selected.}
  \item{bootstrap.ncore}{Number of cores to use for parallelization if the boostrap estimation is selected.}

}

\details{\code{th.quant} is the quantile of the distribution of species density at occurrence sites. 
For example, if \code{th.quant} is set to 0.05, the margin of the distribution is drawn by including 95 percent of the species occurrences, removing the more marginal populations. If \code{th.quant} is set to 0, the margin of the distribution is the minimal envelope including 100 percent of the species occurrences.
The bootstrap estimation of the margin allows representing the uncertainty around this margin. By default it returns a contour covering a 95 percents confidence interval (CI), but you can easily choose a custom CI.}

\value{ a list with $niche.envelope = a raster of the niche envelope including the distribution, $niche.contour = a SpatiaLine object of the margin, $niche.density = a raster of the niche density within the niche envelope.}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

\seealso{\code{\link{ecospat.plot.niche.dyn}}}

\examples{
\donttest{

data(ecospat.testData)

sp1 <- ecospat.testData[ecospat.testData$Bromus_erectus == 1, 1:8] # species occurences
clim <- ecospat.testData[4:8] # environment of the study area

#################################### PCA-ENVIRONMENT ##################################
pca.cal <- ade4::dudi.pca(clim, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

####### projection of the species distribution into the environmental space ######
# predict the scores on the axes
scores.clim <- pca.cal$li[(nrow(occ.sp) + 1):nrow(data), ] # scores for global climate
scores.sp1 <- ade4::suprow(pca.cal, sp1[, 4:8]) # scores for sp1

z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1, R = 100, extend.extent = c(0, 1, 0, 0))

####### estimate the margin ######
z1.margin <- ecospat.margin(z1, th.quant = 0, bootstrap = FALSE) # including all occurrences
z1.95margin <- ecospat.margin(z1, th.quant = 0.05, bootstrap = FALSE) # including 95 percent of the occurrences
z1.bootstrap.margin <- ecospat.margin(z1,
  th.quant = 0,
  bootstrap = TRUE, bootstrap.rep = 100
) # bootstrap estimation of the niche limit

####### plot the margin and its uncertainty ######
plot(z1.margin$niche.density) # plot of the niche density
points(z1.ade$sp) # with species occurences
plot(z1.margin$niche.contour, add = T) # limit of the margin if you include all the distribution
plot(z1.95margin$niche.contour, col = 2, add = T) # limit of the margin if you exclude the 5 percents of the most marginal distribution
plot(z1.bootstrap.margin$niche.contour, col = 2, add = T, lty = "dotted") # limit of the niche based on the 95 percent CI of the bootstrap

plot(z1.bootstrap.margin$niche.envelope) # shows the uncertainty of the niche limit in raster mode
contour(z1.bootstrap.margin$niche.envelope, col = gray(10:1 / 10)) # shows the uncertainty of the niche limit in vector mode

confInt80 <- raster::rasterToPolygons(z1.bootstrap.margin$niche.envelope >= 50, dissolve = TRUE, fun = function(x) {
  x == 1
}) # select a customized confidence interval (here for example 80 percent)
plot(as(confInt80, "SpatialLines"), add = T, lty = "dotted", col = 2)
