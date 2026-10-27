## ----setup, include=FALSE-----------------------------------------------------

local_tmp <- tempdir()

knitr::opts_knit$set(root.dir = local_tmp)

knitr::opts_chunk$set(
  fig.path = file.path(local_tmp, "figures/"),
  cache.path = file.path(local_tmp, "cache/"),
  fig.retina = 1,
  dev = "png",   
  warning = FALSE,
  message = FALSE
)

# Force pandoc to use local directory internally
options(rmarkdown.pandoc.cache = local_tmp)

set.seed(123)

options(rasterTmpDir = local_tmp)

# Load data

## ----load_library-------------------------------------------------------------
library(ecospat)
citation("ecospat")

## ----Test data for the ecospat library----------------------------------------
data(ecospat.testData)
names(ecospat.testData)

## ----Test data for the Niche Overlap Analysis---------------------------------
data(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)
data(ecospat.testNiche.nat)
names(ecospat.testNiche.nat)

## ----Phylogenetic Diversity Analysis------------------------------------------
if(requireNamespace("ape")){
  fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
  tree<-ape::read.tree(fpath)
  tree$tip.label
  plot(tree, cex=0.6)
}

## ----mantel_cor---------------------------------------------------------------
ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, 
                           colvar=3:7, max=1000, nclass=10, nperm=100)

## ----Predictor Variable Selection---------------------------------------------
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)

## ----Number of Predictors-----------------------------------------------------
x <- cor(colvar, method="spearman")
ecospat.npred (x, th=0.75)

## ----Climate Analogy 1--------------------------------------------------------
x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] # A reference dataset

## ----Climate Analogy 2--------------------------------------------------------
ecospat.climan(ref,p)

## ----MESS---------------------------------------------------------------------
x <- ecospat.testData[c(2,3,4:8)]
proj<- x[1:90,] #A projection dataset.
cal<- x[91:300,] #A calibration dataset
mess.object<-ecospat.mess (proj, cal, w="default")
ecospat.plot.mess (mess.object, cex=1, pch=15)

## ----Phylogenetic Diversity---------------------------------------------------
if(requireNamespace("ape")){
  fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
  tree <- ape::read.tree(fpath)
  data <- ecospat.testData[9:52]
  pd<- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = TRUE, average = FALSE, verbose = FALSE )
  plot(pd)
}

## ----data for Niche Quantification--------------------------------------------
library(ade4)
inv <- ecospat.testNiche.inv
nat <- ecospat.testNiche.nat
pca.env <- ade4::dudi.pca(rbind(nat,inv)[,3:10],scannf=F,nf=2) 
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

## ----PCA----------------------------------------------------------------------
# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for the species native distribution
scores.sp.nat <- ade4::suprow(pca.env,nat[which(nat[,11]==1),3:10])$li

# PCA scores for the species invasive distribution
scores.sp.inv <- ade4::suprow(pca.env,inv[which(inv[,11]==1),3:10])$li

# PCA scores for the whole native study area
scores.clim.nat <- ade4::suprow(pca.env,nat[,3:10])$li

# PCA scores for the whole invaded study area
scores.clim.inv <- ade4::suprow(pca.env,inv[,3:10])$li

## ----gridding the native niche------------------------------------------------
# gridding the native niche
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, 
                                       glob1=scores.clim.nat,
                                       sp=scores.sp.nat, R=100,
                                       th.sp=0) 

## ----gridding the invasive niche----------------------------------------------
# gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv,
                                       sp=scores.sp.inv, R=100,
                                       th.sp=0) 

## ----Schoeners D--------------------------------------------------------------
# Compute Schoener's D, index of niche overlap
D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D 
D.overlap

## ----Niche Equivalency Test---------------------------------------------------
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv,rep=10,
                                          intersection = 0.1,
                                          overlap.alternative =  "higher",
                                          expansion.alternative = "lower",
                                          stability.alternative = "higher",
                                          unfilling.alternative = "lower")

## ----Plot Equivalency test----------------------------------------------------
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

## ----Niche Similarity Test----------------------------------------------------
sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=10,
                                          overlap.alternative =  "higher",
                                          expansion.alternative = "lower",
                                          stability.alternative = "higher",
                                          unfilling.alternative = "lower",
                                          intersection = 0.1,
                                          rand.type=1) 

## ----Plot Similarity test-----------------------------------------------------
ecospat.plot.overlap.test(sim.test, "D", "Similarity")

## ----niche.dyn----------------------------------------------------------------
niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv)

## ----Plot niche overlap-------------------------------------------------------
ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")

ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)

## ----Plot Similarity test for niche expansion, stability and unfilling--------
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")

## ----Plot the niche dynamics along one gradient-------------------------------
# gridding the native niche
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]), 
                                         glob1=as.data.frame(nat[,10]),
                                    sp=as.data.frame(nat[which(nat[,11]==1),10]),
                                    R=1000, th.sp=0) 
# gridding the invaded niche
grid.clim.t.inv <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,10]),
                                         glob1=as.data.frame(inv[,10]),
                                         sp=as.data.frame(inv[which(inv[,11]==1),10]), 
                                         R=1000, th.sp=0) 
t.dyn<-ecospat.niche.dyn.index (grid.clim.t.nat, grid.clim.t.inv)
ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0, 
                       interest=2, title= "Niche Overlap", 
                       name.axis1="Average temperature")

## ----co_occ 1-----------------------------------------------------------------
data <- ecospat.testData[c(9:16,54:57)]

## ----co_occ 2-----------------------------------------------------------------
ecospat.co_occurrences (data)

## ----Cscore-------------------------------------------------------------------
data<- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
nperm <- 100
outpath <- getwd()
ecospat.Cscore(data, nperm, outpath)

## ----cor-plot-----------------------------------------------------------------
data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)

## ----Calibration And Evaluation Dataset---------------------------------------
data <- ecospat.testData
caleval <- ecospat.caleval (data = ecospat.testData[53], xy = data[2:3],
                            row.num = 1:nrow(data), nrep = 2, ratio = 0.7,
                            disaggregate = 0.2, pseudoabs = 100, npres = 10,
                            replace = FALSE)
head(caleval)

## ----fit----------------------------------------------------------------------
fit <- ecospat.testData$glm_Saxifraga_oppositifolia

## ----obs----------------------------------------------------------------------
obs<-ecospat.testData$glm_Saxifraga_oppositifolia[which(ecospat.testData$Saxifraga_oppositifolia==1)]


## ----boyce--------------------------------------------------------------------
ecospat.boyce (fit, obs, nclass = 0, window.w = "default", res = 100, 
               PEplot = TRUE)$cor

## ----Accuracy-----------------------------------------------------------------
eval<-ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred<-ecospat.testData[c(73:92)]


## ----ecospat.CommunityEval----------------------------------------------------
CommunityEval<-ecospat.CommunityEval (eval, pred, proba = TRUE, ntir=5,verbose = T)


## ----ESM prep-----------------------------------------------------------------
library(biomod2)

# species
# occurrences
data(ecospat.testNiche.inv)
inv <- ecospat.testNiche.inv
xy <- inv[,1:2]
head(xy)
sp_occ <- inv[11]

# env
current <- inv[3:7]
head(current)

## BIOMOD
t1 <- Sys.time()
sp<-1

## ----ESM format---------------------------------------------------------------
### Formating the data with the BIOMOD_FormatingData() function form the package biomod2

myBiomodData <- biomod2::BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
                                      expl.var = current,
                                      resp.xy = xy,
                                      resp.name = colnames(sp_occ)[sp])

## ----ESM.Modeling, message=FALSE, warning=FALSE-------------------------------
### Calibration of simple bivariate models

# remove insivible(capture.output)) to print output in the console
# this is just to keep the vignette short
invisible(capture.output(
    my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
                                models=c('GLM'),
                                NbRunEval=2,
                                DataSplit=70,
                                weighting.score=c("AUC"),
                                parallel=F)
  ))

## ----ESM.EnsembleModeling-----------------------------------------------------
### Evaluation and average of simple bivariate models to ESMs
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)

## ----ESM.Projection-----------------------------------------------------------
### Projection of simple bivariate models into new space 
my.ESM_proj_current <- ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
                                              new.env=current)

## ----ESM.EnsembleProjection---------------------------------------------------
### Projection of calibrated ESMs into new space 
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
                                                        ESM.EnsembleModeling.output=my.ESM_EF)

## ----spatial pred data--------------------------------------------------------
proba <- ecospat.testData[,73:92]

## ----species richness---------------------------------------------------------
sr <- as.data.frame(rowSums(proba))

## ----SESAM--------------------------------------------------------------------
prr<-ecospat.SESAM.prr(proba, sr)
head(prr)[,1:4]

## ----assemblage pred----------------------------------------------------------
presence<-ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred<-ecospat.testData[c(73:92)]

## ----cons_Score---------------------------------------------------------------
nbpermut <- 100
outpath <- getwd()
ecospat.cons_Cscore(presence, pred, nbpermut, outpath)

