## ----load_library-------------------------------------------------------------
library(ecospat)
citation("ecospat")

## -----------------------------------------------------------------------------
data(ecospat.testData)
names(ecospat.testData)

## -----------------------------------------------------------------------------
data(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)

## -----------------------------------------------------------------------------
data(ecospat.testNiche.nat)
names(ecospat.testNiche.nat)

## -----------------------------------------------------------------------------
fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
fpath
library(ape)
tree<-read.tree(fpath)
tree$tip.label

## ----tree---------------------------------------------------------------------
plot(tree, cex=0.6)

## ----mantel_cor---------------------------------------------------------------
ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, 
                           colvar=3:7, max=1000, nclass=10, nperm=100)

## -----------------------------------------------------------------------------
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)

## -----------------------------------------------------------------------------
x <- cor(colvar, method="spearman")
ecospat.npred (x, th=0.75)

## -----------------------------------------------------------------------------
x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] # A reference dataset

## -----------------------------------------------------------------------------
ecospat.climan(ref,p)

## -----------------------------------------------------------------------------
x <- ecospat.testData[c(2,3,4:8)]
proj<- x[1:90,] #A projection dataset.
cal<- x[91:300,] #A calibration dataset

## -----------------------------------------------------------------------------
mess.object<-ecospat.mess (proj, cal, w="default")

## ----mess---------------------------------------------------------------------
ecospat.plot.mess (mess.object, cex=1, pch=15)

## -----------------------------------------------------------------------------
fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
library(ape)
tree <- read.tree(fpath)
data <- ecospat.testData[9:52]

## -----------------------------------------------------------------------------
pd<- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = TRUE, average = FALSE, verbose = TRUE )

## -----------------------------------------------------------------------------
pd

## ----pd-----------------------------------------------------------------------
plot(pd)

## -----------------------------------------------------------------------------
inv <- ecospat.testNiche.inv

## -----------------------------------------------------------------------------
nat <- ecospat.testNiche.nat

## -----------------------------------------------------------------------------
library(ade4)
pca.env <- dudi.pca(rbind(nat,inv)[,3:10],scannf=F,nf=2) 

## -----------------------------------------------------------------------------
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

## -----------------------------------------------------------------------------
# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,nat[which(nat[,11]==1),3:10])$li

# PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,inv[which(inv[,11]==1),3:10])$li

# PCA scores for the whole native study area
scores.clim.nat <- suprow(pca.env,nat[,3:10])$li

# PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,inv[,3:10])$li

## -----------------------------------------------------------------------------
# gridding the native niche
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, 
                                       glob1=scores.clim.nat,
                                       sp=scores.sp.nat, R=100,
                                       th.sp=0) 

## -----------------------------------------------------------------------------
# gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv,
                                       sp=scores.sp.inv, R=100,
                                       th.sp=0) 

## -----------------------------------------------------------------------------
# Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D 
D.overlap

## -----------------------------------------------------------------------------
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv,rep=10,
                                          intersection = 0.1,
                                          overlap.alternative =  "higher",
                                          expansion.alternative = "lower",
                                          stability.alternative = "higher",
                                          unfilling.alternative = "lower")

## -----------------------------------------------------------------------------
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

## -----------------------------------------------------------------------------
sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=10,
                                          overlap.alternative =  "higher",
                                          expansion.alternative = "lower",
                                          stability.alternative = "higher",
                                          unfilling.alternative = "lower",
                                          intersection = 0.1,
                                          rand.type=1) 

## -----------------------------------------------------------------------------
ecospat.plot.overlap.test(sim.test, "D", "Similarity")

## ----niche.dyn----------------------------------------------------------------
niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = 0.1)

