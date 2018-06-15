pkgname <- "ecospat"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ecospat')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ecospat.CCV.communityEvaluation.bin")
### * ecospat.CCV.communityEvaluation.bin

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.communityEvaluation.bin
### Title: Calculates a range of community evaluation metrics based on
###   different thresholding techniques.
### Aliases: ecospat.CCV.communityEvaluation.bin
### Keywords: ~kwd1 ~kwd2

### ** Examples

#Loading species occurence data and remove empty communities
testData <- ecospat.testData[,c(24,34,43,45,48,53,55:58,60:63,65:66,68:71)]
sp.data <- testData[which(rowSums(testData)>0), sort(colnames(testData))]

#Loading environmental data
env.data <- ecospat.testData[which(rowSums(testData)>0),4:8]

#Coordinates for all sites
xy <- ecospat.testData[which(rowSums(testData)>0),2:3]

#Running all the models for all species
myCCV.Models <- ecospat.CCV.modeling(sp.data = sp.data,
                                     env.data = env.data,
                                     xy = xy,
                                     NbRunEval = 5,
                                     minNbPredictors = 10,
                                     VarImport = 3)
                                     
#Thresholding all the predictions and calculating the community evaluation metrics
myCCV.communityEvaluation.bin <- ecospat.CCV.communityEvaluation.bin(
      ccv.modeling.data = myCCV.Models, 
      thresholds = c('MAX.KAPPA', 'MAX.ROC','PS_SDM'),
      community.metrics= c('SR.deviation','Sorensen'),
      parallel = FALSE,
      cpus = 4)



cleanEx()
nameEx("ecospat.CCV.communityEvaluation.prob")
### * ecospat.CCV.communityEvaluation.prob

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.communityEvaluation.prob
### Title: Evaluates community predictions directly on the probabilities
###   (i.e., threshold independent)
### Aliases: ecospat.CCV.communityEvaluation.prob
### Keywords: ~kwd1 ~kwd2

### ** Examples

#Loading species occurence data and remove empty communities
testData <- ecospat.testData[,c(24,34,43,45,48,53,55:58,60:63,65:66,68:71)]
sp.data <- testData[which(rowSums(testData)>0), sort(colnames(testData))]

#Loading environmental data
env.data <- ecospat.testData[which(rowSums(testData)>0),4:8]

#Coordinates for all sites
xy <- ecospat.testData[which(rowSums(testData)>0),2:3]

#Running all the models for all species
myCCV.Models <- ecospat.CCV.modeling(sp.data = sp.data,
                                     env.data = env.data,
                                     xy = xy,
                                     NbRunEval = 5,
                                     minNbPredictors = 10,
                                     VarImport = 3)
                                     
#Calculating the probabilistic community metrics
myCCV.communityEvaluation.prob <- ecospat.CCV.communityEvaluation.prob(
      ccv.modeling.data = myCCV.Models,
      community.metrics = c('SR.deviation','community.AUC','probabilistic.Sorensen'),
      se.th = 0.02, 
      parallel = FALSE,
      cpus = 4)



cleanEx()
nameEx("ecospat.CCV.createDataSplitTable")
### * ecospat.CCV.createDataSplitTable

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.createDataSplitTable
### Title: Creates a DataSplitTable for usage in ecospat.ccv.modeling.
### Aliases: ecospat.CCV.createDataSplitTable
### Keywords: ~kwd1 ~kwd2

### ** Examples

#Creating a DataSplitTable for 200 sites, 25 runs with an 
#80/20 calibration/evaluation cross-validation

DataSplitTable <- ecospat.CCV.createDataSplitTable(NbSites = 200, 
                                                   NbRunEval=25, 
                                                   DataSplit=80, 
                                                   validation.method='cross-validation')
                                                   
#Loading species occurence data and remove empty communities
testData <- ecospat.testData[,c(24,34,43,45,48,53,55:58,60:63,65:66,68:71)]
sp.data <- testData[which(rowSums(testData)>0), sort(colnames(testData))]

#Creating a DataSplitTable based on species data directly
DataSplitTable <- ecospat.CCV.createDataSplitTable(NbRunEval = 20,
                                                   DataSplit = 70,
                                                   validation.method = "cross-validation",
                                                   NbSites = NULL,
                                                   sp.data = sp.data, 
                                                   minNbPresence = 15, 
                                                   minNbAbsences = 15, 
                                                   maxNbTry = 250)



cleanEx()
nameEx("ecospat.CCV.modeling")
### * ecospat.CCV.modeling

flush(stderr()); flush(stdout())

### Name: ecospat.CCV.modeling
### Title: Runs indivudual species distribuion models with SDMs or ESMs
### Aliases: ecospat.CCV.modeling
### Keywords: ~kwd1 ~kwd2

### ** Examples

#Loading species occurence data and remove empty communities
testData <- ecospat.testData[,c(24,34,43,45,48,53,55:58,60:63,65:66,68:71)]
sp.data <- testData[which(rowSums(testData)>0), sort(colnames(testData))]

#Loading environmental data
env.data <- ecospat.testData[which(rowSums(testData)>0),4:8]

#Coordinates for all sites
xy <- ecospat.testData[which(rowSums(testData)>0),2:3]

#Running all the models for all species
myCCV.Models <- ecospat.CCV.modeling(sp.data = sp.data,
                                     env.data = env.data,
                                     xy = xy,
                                     NbRunEval = 5,
                                     minNbPredictors = 10,
                                     VarImport = 3)



cleanEx()
nameEx("ecospat.CommunityEval")
### * ecospat.CommunityEval

flush(stderr()); flush(stdout())

### Name: ecospat.CommunityEval
### Title: Community Evaluation
### Aliases: ecospat.CommunityEval

### ** Examples

## Not run: 
##D eval <- Data[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
##D pred <- Data[c(73:92)]
##D 
##D ecospat.CommunityEval (eval, pred, proba=TRUE, ntir=10)
## End(Not run)


cleanEx()
nameEx("ecospat.Cscore")
### * ecospat.Cscore

flush(stderr()); flush(stdout())

### Name: ecospat.Cscore
### Title: Pairwise co-occurrence Analysis with calculation of the C-score
###   index.
### Aliases: ecospat.Cscore

### ** Examples

## Not run: 
##D data<- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
##D nperm <- 10000
##D outpath <- getwd()
##D ecospat.Cscore(data, nperm, outpath)
##D 
## End(Not run)


cleanEx()
nameEx("ecospat.ESM.EnsembleModeling")
### * ecospat.ESM.EnsembleModeling

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.EnsembleModeling
### Title: Ensamble of Small Models: Evaluates and Averages Simple
###   Bivariate Models To ESMs
### Aliases: ecospat.ESM.EnsembleModeling

### ** Examples

   ## Not run: 
##D # Loading test data
##D inv <- ecospat.testNiche.inv
##D 
##D # species occurrences
##D xy <- inv[,1:2]
##D sp_occ <- inv[11]
##D 
##D # env
##D current <- inv[3:10]
##D 
##D 
##D 
##D ### Formating the data with the BIOMOD_FormatingData() function from the package biomod2
##D 
##D sp <- 1
##D myBiomodData <- BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
##D                                       expl.var = current,
##D                                       resp.xy = xy,
##D                                       resp.name = colnames(sp_occ)[sp])
##D 
##D 
##D ### Calibration of simple bivariate models
##D my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
##D                                 models=c('GLM','RF'),
##D                                 NbRunEval=2,
##D                                 DataSplit=70,
##D                                 weighting.score=c("AUC"),
##D                                 parallel=FALSE)  
##D 
##D 
##D ### Evaluation and average of simple bivariate models to ESMs
##D my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)
##D 
##D ### Projection of simple bivariate models into new space 
##D my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
##D                                             new.env=current)
##D 
##D ### Projection of calibrated ESMs into new space 
##D my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
##D                                                         ESM.EnsembleModeling.output=my.ESM_EF)
##D 
##D ## get the model performance of ESMs 
##D my.ESM_EF$ESM.evaluations
##D ## get the weights of the single bivariate models used to build the ESMs
##D my.ESM_EF$weights
## End(Not run)


cleanEx()
nameEx("ecospat.ESM.EnsembleProjection")
### * ecospat.ESM.EnsembleProjection

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.EnsembleProjection
### Title: Ensamble of Small Models: Projects Calibrated ESMs Into New
###   Space Or Time.
### Aliases: ecospat.ESM.EnsembleProjection

### ** Examples

   ## Not run: 
##D # Loading test data for the niche dynamics analysis in the invaded range
##D inv <- ecospat.testNiche.inv
##D 
##D # species occurrences
##D xy <- inv[,1:2]
##D sp_occ <- inv[11]
##D 
##D # env
##D current <- inv[3:10]
##D 
##D 
##D 
##D ### Formating the data with the BIOMOD_FormatingData() function form the package biomod2
##D setwd(path.wd)
##D t1 <- Sys.time()
##D sp <- 1
##D myBiomodData <- BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
##D                                       expl.var = current,
##D                                       resp.xy = xy,
##D                                       resp.name = colnames(sp_occ)[sp])
##D 
##D myBiomodOption <- Print_Default_ModelingOptions()
##D 
##D 
##D ### Calibration of simple bivariate models
##D my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
##D                                 models=c('GLM','RF'),
##D                                 models.options=myBiomodOption,
##D                                 NbRunEval=2,
##D                                 DataSplit=70,
##D                                 weighting.score=c("AUC"),
##D                                 parallel=FALSE)  
##D 
##D 
##D ### Evaluation and average of simple bivariate models to ESMs
##D my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)
##D 
##D ### Projection of simple bivariate models into new space 
##D my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
##D                                             new.env=current)
##D 
##D ### Projection of calibrated ESMs into new space 
##D my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
##D                                                         ESM.EnsembleModeling.output=my.ESM_EF)
##D 
##D ## get the model performance of ESMs 
##D my.ESM_EF$ESM.evaluations
##D ## get the weights of the single bivariate models used to build the ESMs
##D my.ESM_EF$weights
## End(Not run)


cleanEx()
nameEx("ecospat.ESM.Modeling")
### * ecospat.ESM.Modeling

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.Modeling
### Title: Ensamble of Small Models: Calibration of Simple Bivariate Models
### Aliases: ecospat.ESM.Modeling

### ** Examples

   ## Not run: 
##D # Loading test data
##D inv <- ecospat.testNiche.inv
##D 
##D # species occurrences
##D xy <- inv[,1:2]
##D sp_occ <- inv[11]
##D 
##D # env
##D current <- inv[3:10]
##D 
##D 
##D 
##D ### Formating the data with the BIOMOD_FormatingData() function from the package biomod2
##D sp <- 1
##D myBiomodData <- BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
##D                                       expl.var = current,
##D                                       resp.xy = xy,
##D                                       resp.name = colnames(sp_occ)[sp])
##D 
##D ### Calibration of simple bivariate models
##D my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
##D                                 models=c('GLM','RF'),
##D                                 NbRunEval=2,
##D                                 DataSplit=70,
##D                                 Prevalence=0.5
##D                                 weighting.score=c("AUC"),
##D                                 parallel=FALSE)  
##D 
##D 
##D ### Evaluation and average of simple bivariate models to ESMs
##D my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)
##D 
##D ### Projection of simple bivariate models into new space 
##D my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
##D                                             new.env=current)
##D 
##D ### Projection of calibrated ESMs into new space 
##D my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
##D                                                         ESM.EnsembleModeling.output=my.ESM_EF)
##D 
##D ## get the model performance of ESMs 
##D my.ESM_EF$ESM.evaluations
##D ## get the weights of the single bivariate models used to build the ESMs
##D my.ESM_EF$weights
## End(Not run)


cleanEx()
nameEx("ecospat.ESM.Projection")
### * ecospat.ESM.Projection

flush(stderr()); flush(stdout())

### Name: ecospat.ESM.Projection
### Title: Ensamble of Small Models: Projects Simple Bivariate Models Into
###   New Space Or Time
### Aliases: ecospat.ESM.Projection

### ** Examples

   ## Not run: 
##D # Loading test data
##D inv <- ecospat.testNiche.inv
##D 
##D # species occurrences
##D xy <- inv[,1:2]
##D sp_occ <- inv[11]
##D 
##D # env
##D current <- inv[3:10]
##D 
##D 
##D 
##D ### Formating the data with the BIOMOD_FormatingData() function from the package biomod2
##D        
##D sp <- 1
##D myBiomodData <- BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
##D                                       expl.var = current,
##D                                       resp.xy = xy,
##D                                       resp.name = colnames(sp_occ)[sp])
##D 
##D ### Calibration of simple bivariate models
##D my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
##D                                 models=c('GLM','RF'),
##D                                 NbRunEval=2,
##D                                 DataSplit=70,
##D                                 weighting.score=c("AUC"),
##D                                 parallel=FALSE)  
##D 
##D 
##D ### Evaluation and average of simple bivariate models to ESMs
##D my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)
##D 
##D ### Projection of simple bivariate models into new space 
##D my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
##D                                             new.env=current)
##D 
##D ### Projection of calibrated ESMs into new space 
##D my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
##D                                                         ESM.EnsembleModeling.output=my.ESM_EF)
##D 
##D ## get the model performance of ESMs 
##D my.ESM_EF$ESM.evaluations
##D ## get the weights of the single bivariate models used to build the ESMs
##D my.ESM_EF$weights
## End(Not run)


cleanEx()
nameEx("ecospat.Epred")
### * ecospat.Epred

flush(stderr()); flush(stdout())

### Name: ecospat.Epred
### Title: Prediction Mean
### Aliases: ecospat.Epred

### ** Examples

x <- ecospat.testData[c(92,96)]
mean <- ecospat.Epred (x, w=rep(1,ncol(x)), th=0.5)



cleanEx()
nameEx("ecospat.SESAM.prr")
### * ecospat.SESAM.prr

flush(stderr()); flush(stdout())

### Name: ecospat.SESAM.prr
### Title: SESAM Probability Ranking Rule
### Aliases: ecospat.SESAM.prr

### ** Examples

proba <- ecospat.testData[,73:92]
sr <- as.data.frame(rowSums(proba))
ecospat.SESAM.prr(proba, sr)




cleanEx()
nameEx("ecospat.adj.D2")
### * ecospat.adj.D2

flush(stderr()); flush(stdout())

### Name: ecospat.adj.D2.glm
### Title: Calculate An Adjusted D2
### Aliases: ecospat.adj.D2.glm

### ** Examples


glm.obj<-glm(Achillea_millefolium~ddeg+mind+srad+slp+topo, 
family = binomial, data=ecospat.testData)

ecospat.adj.D2.glm(glm.obj)




cleanEx()
nameEx("ecospat.binary.model")
### * ecospat.binary.model

flush(stderr()); flush(stdout())

### Name: ecospat.binary.model
### Title: Generate Binary Models
### Aliases: ecospat.binary.model

### ** Examples

library(dismo)

# get predictor variables
fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
                     pattern='grd', full.names=TRUE )
predictors <- stack(fnames)


# file with presence points
occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
colnames(occ) <- c("x","y")

# fit a domain model, biome is a categorical variable

do <- domain(predictors, occ, factors='biome')

# predict to entire dataset
pred <- predict(do, predictors) 

plot(pred)
points(occ)


# use MPA to convert suitability to binary map (90% of occurrences encompass by binary map)
mpa.cutoff <- ecospat.mpa(pred,occ)

pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)

plot(pred.bin.mpa)
points(occ)



cleanEx()
nameEx("ecospat.boyce")
### * ecospat.boyce

flush(stderr()); flush(stdout())

### Name: ecospat.boyce
### Title: Calculate Boyce Index
### Aliases: ecospat.boyce

### ** Examples

obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.boyce (fit = ecospat.testData$glm_Saxifraga_oppositifolia , obs, nclass=0, 
window.w="default", res=100, PEplot = TRUE)



cleanEx()
nameEx("ecospat.calculate.pd")
### * ecospat.calculate.pd

flush(stderr()); flush(stdout())

### Name: ecospat.calculate.pd
### Title: Calculate Phylogenetic Diversity Measures
### Aliases: ecospat.calculate.pd

### ** Examples

fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
tree <-read.tree(fpath)
data <- ecospat.testData[9:52] 

pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = FALSE, 
average = FALSE, verbose = TRUE )

plot(pd)



cleanEx()
nameEx("ecospat.caleval")
### * ecospat.caleval

flush(stderr()); flush(stdout())

### Name: ecospat.caleval
### Title: Calibration And Evaluation Dataset
### Aliases: ecospat.caleval

### ** Examples

data <- ecospat.testData
caleval <- ecospat.caleval (data = ecospat.testData[53], xy = data[2:3], row.num = 1:nrow(data), 
nrep = 2, ratio = 0.7, disaggregate = 0.2, pseudoabs = 100, npres = 10, replace = FALSE)
caleval



cleanEx()
nameEx("ecospat.climan")
### * ecospat.climan

flush(stderr()); flush(stdout())

### Name: ecospat.climan
### Title: A climate analogy setection tool for the modeling of species
###   distributions
### Aliases: ecospat.climan

### ** Examples

x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] #A reference dataset
ecospat.climan(ref,p)




cleanEx()
nameEx("ecospat.co_occurrences")
### * ecospat.co_occurrences

flush(stderr()); flush(stdout())

### Name: ecospat.co_occurrences
### Title: Species Co-Occurrences
### Aliases: ecospat.co_occurrences

### ** Examples

## Not run: 
##D matrix <- ecospat.testData[c(9:16,54:57)]
##D ecospat.co_occurrences (data=matrix)
## End(Not run)


cleanEx()
nameEx("ecospat.cohen.kappa")
### * ecospat.cohen.kappa

flush(stderr()); flush(stdout())

### Name: ecospat.cohen.kappa
### Title: Cohen's Kappa
### Aliases: ecospat.cohen.kappa

### ** Examples

Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
th <- 0.39 # threshold
xtab <- table(Pred >= th, Sp.occ)

ecospat.cohen.kappa(xtab)



cleanEx()
nameEx("ecospat.cons_Cscore")
### * ecospat.cons_Cscore

flush(stderr()); flush(stdout())

### Name: ecospat.cons_Cscore
### Title: Constrained Co-Occurrence Analysis.
### Aliases: ecospat.cons_Cscore

### ** Examples

## Not run: 
##D presence <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
##D pred <- ecospat.testData[c(73:92)]
##D nperm <- 10000
##D outpath <- getwd()
##D ecospat.cons_Cscore(presence, pred, nperm, outpath)
## End(Not run)


cleanEx()
nameEx("ecospat.cor.plot")
### * ecospat.cor.plot

flush(stderr()); flush(stdout())

### Name: ecospat.cor.plot
### Title: Correlation Plot
### Aliases: ecospat.cor.plot

### ** Examples

data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)



cleanEx()
nameEx("ecospat.cv.example")
### * ecospat.cv.example

flush(stderr()); flush(stdout())

### Name: ecospat.cv.example
### Title: Cross Validation Example Function
### Aliases: ecospat.cv.example

### ** Examples

## Not run: 
##D  
##D ecospat.cv.example ()
## End(Not run)



cleanEx()
nameEx("ecospat.cv.gbm")
### * ecospat.cv.gbm

flush(stderr()); flush(stdout())

### Name: ecospat.cv.gbm
### Title: GBM Cross Validation
### Aliases: ecospat.cv.gbm

### ** Examples

## Not run: 
##D gbm <- ecospat.cv.gbm (gbm.obj= get ("gbm.Agrostis_capillaris", envir=ecospat.env), 
##D ecospat.testData, K=10, cv.lim=10, jack.knife=FALSE)
## End(Not run)



cleanEx()
nameEx("ecospat.cv.glm")
### * ecospat.cv.glm

flush(stderr()); flush(stdout())

### Name: ecospat.cv.glm
### Title: GLM Cross Validation
### Aliases: ecospat.cv.glm

### ** Examples

## Not run: 
##D glm <- ecospat.cv.glm (glm.obj = get ("glm.Agrostis_capillaris", envir=ecospat.env), 
##D K=10, cv.lim=10, jack.knife=FALSE)
## End(Not run)



cleanEx()
nameEx("ecospat.cv.me")
### * ecospat.cv.me

flush(stderr()); flush(stdout())

### Name: ecospat.cv.me
### Title: Maxent Cross Validation
### Aliases: ecospat.cv.me

### ** Examples


## Not run: 
##D me <- ecospat.cv.me(ecospat.testData, names(ecospat.testData)[53], 
##D names(ecospat.testData)[4:8], K = 10, cv.lim = 10, jack.knife = FALSE)
## End(Not run)



cleanEx()
nameEx("ecospat.cv.rf")
### * ecospat.cv.rf

flush(stderr()); flush(stdout())

### Name: ecospat.cv.rf
### Title: RandomForest Cross Validation
### Aliases: ecospat.cv.rf

### ** Examples

## Not run: 
##D rf <- ecospat.cv.rf(get("rf.Agrostis_capillaris", envir = ecospat.env), 
##D ecospat.testData[, c(53, 4:8)], K = 10, cv.lim = 10, jack.knife = FALSE)
## End(Not run)



cleanEx()
nameEx("ecospat.env")
### * ecospat.env

flush(stderr()); flush(stdout())

### Name: ecospat.env
### Title: Package Environment
### Aliases: ecospat.env

### ** Examples

ls(envir=ecospat.env)



cleanEx()
nameEx("ecospat.grid.clim.dyn")
### * ecospat.grid.clim.dyn

flush(stderr()); flush(stdout())

### Name: ecospat.grid.clim.dyn
### Title: Dynamic Occurrence Densities Grid
### Aliases: ecospat.grid.clim.dyn

### ** Examples

## Not run: 
##D spp <- ecospat.testNiche
##D clim <- ecospat.testData[2:8]
##D 
##D occ.sp_test <- na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=2:3,colspkept=1:3,dfvar=clim,
##D colvarxy=1:2,colvar="all",resolution=25))
##D 
##D occ.sp<-cbind(occ.sp_test,spp[,4]) #add species names
##D 
##D # list of species
##D sp.list<-levels(occ.sp[,1])
##D sp.nbocc<-c()
##D 
##D for (i in 1:length(sp.list)){sp.nbocc<-c(sp.nbocc,length(which(occ.sp[,1] == sp.list[i])))} 
##D #calculate the nb of occurences per species
##D 
##D sp.list <- sp.list[sp.nbocc>4] # remove species with less than 5 occurences
##D nb.sp <- length(sp.list) #nb of species
##D ls()
##D # selection of variables to include in the analyses 
##D # try with all and then try only worldclim Variables
##D Xvar <- c(3:7)
##D nvar <- length(Xvar)
##D 
##D #number of interation for the tests of equivalency and similarity
##D iterations <- 100
##D #resolution of the gridding of the climate space
##D R <- 100
##D #################################### PCA-ENVIRONMENT ##################################
##D data<-rbind(occ.sp[,Xvar+1],clim[,Xvar]) 
##D w <- c(rep(0,nrow(occ.sp)),rep(1,nrow(clim)))
##D pca.cal <- dudi.pca(data, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
##D 
##D ####### selection of species ######
##D sp.list
##D sp.combn <- combn(1:2,2)
##D 
##D for(i in 1:ncol(sp.combn)) {
##D   row.sp1 <- which(occ.sp[,1] == sp.list[sp.combn[1,i]]) # rows in data corresponding to sp1
##D   row.sp2 <- which(occ.sp[,1] == sp.list[sp.combn[2,i]]) # rows in data corresponding to sp2
##D   name.sp1 <- sp.list[sp.combn[1,i]]
##D   name.sp2 <- sp.list[sp.combn[2,i]]
##D   # predict the scores on the axes
##D   scores.clim <- pca.cal$li[(nrow(occ.sp)+1):nrow(data),]  #scores for global climate
##D   scores.sp1 <- pca.cal$li[row.sp1,]					#scores for sp1
##D   scores.sp2 <- pca.cal$li[row.sp2,]					#scores for sp2
##D }
##D # calculation of occurence density and test of niche equivalency and similarity 
##D z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
##D z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2,R=100)
## End(Not run)


cleanEx()
nameEx("ecospat.makeDataFrame")
### * ecospat.makeDataFrame

flush(stderr()); flush(stdout())

### Name: ecospat.makeDataFrame
### Title: Make Data Frame
### Aliases: ecospat.makeDataFrame

### ** Examples

## Not run: 
##D files <- list.files(path=paste(system.file(package="dismo"),
##D                                '/ex', sep=''), pattern='grd', full.names=TRUE )
##D predictors <- raster::stack(files[c(9,1:8)])   #file 9 has more NA values than
##D # the other files, this is why we choose it as the first layer (see ?randomPoints)
##D 
##D file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
##D bradypus <- read.table(file, header=TRUE, sep=',')[,c(2,3,1)]
##D head(bradypus)
##D 
##D random.spec <- cbind(as.data.frame(randomPoints(predictors,50,extf=1)),species="randomSpec")
##D colnames(random.spec)[1:2] <- c("lon","lat")
##D 
##D spec.list <- rbind(bradypus, random.spec)
##D 
##D df <- ecospat.makeDataFrame(spec.list, expl.var=predictors, n=5000)
##D head(df)
##D 
##D plot(predictors[[1]])
##D points(df[df$Bradypus.variegatus==1, c('x','y')])
##D points(df[df$randomSpec==1, c('x','y')], col="red")
##D 
## End(Not run)



cleanEx()
nameEx("ecospat.mantel.correlogram")
### * ecospat.mantel.correlogram

flush(stderr()); flush(stdout())

### Name: ecospat.mantel.correlogram
### Title: Mantel Correlogram
### Aliases: ecospat.mantel.correlogram

### ** Examples

ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, colvar=3:7, 
max=1000, nclass=10, nperm=100)



cleanEx()
nameEx("ecospat.max.kappa")
### * ecospat.max.kappa

flush(stderr()); flush(stdout())

### Name: ecospat.max.kappa
### Title: Maximum Kappa
### Aliases: ecospat.max.kappa

### ** Examples


   ## Not run: 
##D Pred <- ecospat.testData$glm_Agrostis_capillaris
##D Sp.occ <- ecospat.testData$Agrostis_capillaris
##D kappa100 <- ecospat.max.kappa(Pred, Sp.occ)
##D    
## End(Not run)




cleanEx()
nameEx("ecospat.max.tss")
### * ecospat.max.tss

flush(stderr()); flush(stdout())

### Name: ecospat.max.tss
### Title: Maximum TSS
### Aliases: ecospat.max.tss
### Keywords: file

### ** Examples


Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
TSS100 <- ecospat.max.tss(Pred, Sp.occ)



cleanEx()
nameEx("ecospat.maxentvarimport")
### * ecospat.maxentvarimport

flush(stderr()); flush(stdout())

### Name: ecospat.maxentvarimport
### Title: Maxent Variable Importance
### Aliases: ecospat.maxentvarimport

### ** Examples

## Not run: 
##D model <- get ("me.Achillea_millefolium", envir=ecospat.env)
##D dfvar <- ecospat.testData[4:8]
##D nperm <- 5
##D ecospat.maxentvarimport (model, cal, nperm)
## End(Not run)



cleanEx()
nameEx("ecospat.mdr")
### * ecospat.mdr

flush(stderr()); flush(stdout())

### Name: ecospat.mdr
### Title: Minimum Dispersal Routes)
### Aliases: ecospat.mdr

### ** Examples

   ## Not run: 
##D library(maps)
##D data<- ecospat.testMdr
##D 
##D fixed.sources.rows<-order(data$date)[1:2] #first introductions
##D 
##D #plot observed situation
##D plot(data[,2:1],pch=15,cex=0.5)
##D points(data[fixed.sources.rows,2:1],pch=19,col="red")
##D text(data[,2]+0.5,data[,1]+0.5,data[,3],cex=0.5)
##D map(add=T)
##D 
##D # mca 
##D obs<-ecospat.mdr(data=data,
##D xcol=2,
##D ycol=1,
##D datecol=3,
##D mode="min",
##D rep=100,
##D mean.date.error=1,
##D fixed.sources.rows)
##D 
##D #plot results
##D lwd<-(obs[[1]]$bootstrap.value)
##D x11();plot(obs[[1]][,3:4],type="n",xlab="longitude",ylab="latitude")
##D arrows(obs[[1]][,1],obs[[1]][,2],obs[[1]][,3],obs[[1]][,4],length = 0.05,lwd=lwd*2)
##D map(add=T)
##D points(data[fixed.sources.rows,2:1],pch=19,col="red")
##D text(data[fixed.sources.rows,2]+0.5,data[fixed.sources.rows,1]+0.5,data[fixed.sources.rows,3],
##D cex=1,col="red")
##D title(paste("total routes length : ",
##D round(obs[[2]],2)," Deg","\n","median dispersal rate : ",
##D round(obs[[3]],2)," Deg/year","\n","number of outcoming nodes : ",
##D obs[[4]]))
## End(Not run)



cleanEx()
nameEx("ecospat.mess")
### * ecospat.mess

flush(stderr()); flush(stdout())

### Name: ecospat.mess
### Title: MESS
### Aliases: ecospat.mess

### ** Examples

x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset.
cal <- x[91:300,] #A calibration dataset

#Create a MESS object 
mess.object <- ecospat.mess (proj, cal, w="default")

#Plot MESS 
ecospat.plot.mess (mess.object, cex=1, pch=15)



cleanEx()
nameEx("ecospat.meva.table")
### * ecospat.meva.table

flush(stderr()); flush(stdout())

### Name: ecospat.meva.table
### Title: Model Evaluation For A Given Threshold Value
### Aliases: ecospat.meva.table
### Keywords: file

### ** Examples


Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris

meva <- ecospat.meva.table (Pred, Sp.occ, 0.39)



cleanEx()
nameEx("ecospat.migclim")
### * ecospat.migclim

flush(stderr()); flush(stdout())

### Name: ecospat.migclim
### Title: Implementing Dispersal Into Species Distribution Models
### Aliases: ecospat.migclim

### ** Examples

## Not run: 
##D ecospat.migclim()
##D ### Some example data files can be downloaded from the following web page:
##D ### http://www.unil.ch/ecospat/page89413.html
##D ###
##D ### Run the example as follows (set the current working directory to the
##D ### folder where the example data files are located):
##D ###
##D data(MigClim.testData)
##D ### Run MigClim with a data frame type input.
##D n <- MigClim.migrate (iniDist=MigClim.testData[,1:3],
##D hsMap=MigClim.testData[,4:8], rcThreshold=500,
##D envChgSteps=5, dispSteps=5, dispKernel=c(1.0,0.4,0.16,0.06,0.03),
##D barrier=MigClim.testData[,9], barrierType="strong",
##D iniMatAge=1, propaguleProd=c(0.01,0.08,0.5,0.92),
##D lddFreq=0.1, lddMinDist=6, lddMaxDist=15,
##D simulName="MigClimTest", replicateNb=1, overWrite=TRUE,
##D testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)
## End(Not run)




cleanEx()
nameEx("ecospat.mpa")
### * ecospat.mpa

flush(stderr()); flush(stdout())

### Name: ecospat.mpa
### Title: Minimal Predicted Area
### Aliases: ecospat.mpa

### ** Examples

obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.mpa(obs)
ecospat.mpa(obs,perc=1) ## 100 percent of the presences encompassed



cleanEx()
nameEx("ecospat.niche.dynIndexProjGeo")
### * ecospat.niche.dynIndexProjGeo

flush(stderr()); flush(stdout())

### Name: ecospat.niche.dynIndexProjGeo
### Title: Projection of niche dynamic indices to the Geography
### Aliases: ecospat.niche.dynIndexProjGeo

### ** Examples

## Not run: 
##D 
##D library(raster)
##D 
##D spp <- ecospat.testNiche
##D xy.sp1<-subset(spp,species=="sp1")[2:3] #Bromus_erectus
##D xy.sp2<-subset(spp,species=="sp3")[2:3] #Daucus_carota
##D 
##D ?ecospat.testEnvRaster
##D load(system.file("extdata", "ecospat.testEnvRaster.Rdata", package="ecospat"))
##D 
##D env.sp1<-extract(env,xy.sp1)
##D env.sp2<-extract(env,xy.sp2)
##D env.bkg<-na.exclude(values(env))
##D 
##D #################################### PCA-ENVIRONMENT ##################################
##D 
##D pca.cal <- dudi.pca(env.bkg, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
##D 
##D # predict the scores on the axes
##D scores.bkg <- pca.cal$li  #scores for background climate
##D scores.sp1 <- suprow(pca.cal,env.sp1)$lisup				#scores for sp1
##D scores.sp2 <- suprow(pca.cal,env.sp2)$lisup				#scores for sp2
##D 
##D # calculation of occurence density (niche z)
##D   
##D z1 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg, scores.sp1,R=100)
##D z2 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg, scores.sp2,R=100)
##D 
##D plot(z1$z.uncor)
##D points(scores.sp1)
##D 
##D plot(z2$z.uncor)
##D points(scores.sp2)
##D 
##D ecospat.niche.overlap(z1,z2 ,cor=T)
##D 
##D #################################### stability S in space ##################################
##D 
##D geozS<-ecospat.niche.dynIndexProjGeo(z1,z2,env,index="stability")
##D 
##D plot(geozS,main="Stability")
##D points(xy.sp1,col="red")
##D points(xy.sp2,col="blue")
## End(Not run)


cleanEx()
nameEx("ecospat.niche.zProjGeo")
### * ecospat.niche.zProjGeo

flush(stderr()); flush(stdout())

### Name: ecospat.niche.zProjGeo
### Title: Projection of Occurrence Densities to the Geography
### Aliases: ecospat.niche.zProjGeo

### ** Examples

## Not run: 
##D 
##D library(raster)
##D 
##D spp <- ecospat.testNiche
##D xy.sp1<-subset(spp,species=="sp1")[2:3] #Bromus_erectus
##D 
##D load(system.file("extdata", "ecospat.testEnvRaster.Rdata", package="ecospat"))
##D #?ecospat.testEnvRaster
##D 
##D env.sp1<-extract(env,xy.sp1)
##D env.bkg<-na.exclude(values(env))
##D 
##D #################################### PCA-ENVIRONMENT ##################################
##D 
##D pca.cal <- dudi.pca(env.bkg, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
##D 
##D # predict the scores on the axes
##D scores.bkg <- pca.cal$li  #scores for background climate
##D scores.sp1 <- suprow(pca.cal,env.sp1)$lisup				#scores for sp1
##D 
##D # calculation of occurence density (niche z)
##D   
##D z1 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg, scores.sp1,R=100)
##D 
##D plot(z1$z.uncor)
##D points(scores.sp1)
##D 
##D #################################### occurrence density in space ##################################
##D 
##D # sp1
##D geoz1<-ecospat.niche.zProjGeo(z1,env)
##D plot(geoz1,main="z1")
##D points(xy.sp1)
## End(Not run)


cleanEx()
nameEx("ecospat.npred")
### * ecospat.npred

flush(stderr()); flush(stdout())

### Name: ecospat.npred
### Title: Number Of Predictors
### Aliases: ecospat.npred

### ** Examples

colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)



cleanEx()
nameEx("ecospat.occ.desaggregation")
### * ecospat.occ.desaggregation

flush(stderr()); flush(stdout())

### Name: ecospat.occ.desaggregation
### Title: Species Occurrences Desaggregation
### Aliases: ecospat.occ.desaggregation

### ** Examples


## Not run: 
##D spp <- ecospat.testNiche
##D colnames(spp)[2:3] <- c('x','y')
##D sp1 <- spp[1:32,2:3]
##D 
##D occ.sp1 <- ecospat.occ.desaggregation(xy=sp1, min.dist=500, by=NULL)
##D occ.all.sp <- ecospat.occ.desaggregation(xy=spp, min.dist=500, by='Spp')
## End(Not run)



cleanEx()
nameEx("ecospat.occupied.patch")
### * ecospat.occupied.patch

flush(stderr()); flush(stdout())

### Name: ecospat.occupied.patch
### Title: Extract occupied patches of a species in geographic space.)
### Aliases: ecospat.occupied.patch
### Keywords: file

### ** Examples


## Not run: 
##D 
##D library(dismo)
##D 
##D 
##D library(dismo)
##D 
##D 
##D # only run if the maxent.jar file is available, in the right folder
##D jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
##D 
##D # checking if maxent can be run (normally not part of your script)
##D file.exists(jar)
##D require(rJava))
##D 
##D # get predictor variables
##D fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
##D                      pattern='grd', full.names=TRUE )
##D predictors <- stack(fnames)
##D #plot(predictors)
##D 
##D # file with presence points
##D occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
##D occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
##D colnames(occ) <- c("x","y")
##D occ <- ecospat.occ.desaggregation(occ,min.dist=1)
##D 
##D # fit model, biome is a categorical variable
##D me <- maxent(predictors, occ, factors='biome')
##D 
##D # predict to entire dataset
##D pred <- predict(me, predictors) 
##D 
##D plot(pred)
##D points(occ)
##D 
##D 
##D # use MPA to convert suitability to binary map
##D mpa.cutoff <- ecospat.mpa(pred,occ)
##D 
##D pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)
##D names(pred.bin.mpa) <- "me.mpa"
##D pred.bin.arbitrary <- ecospat.binary.model(pred,0.5)
##D names(pred.bin.arbitrary) <- "me.arbitrary"
##D 
##D 
##D mpa.ocp  <- ecospat.occupied.patch(pred.bin.mpa,occ)
##D arbitrary.ocp  <- ecospat.occupied.patch(pred.bin.arbitrary,occ)
##D 
##D par(mfrow=c(1,2))
##D plot(mpa.ocp) ## occupied patches: green area
##D points(occ,col="red",cex=0.5,pch=19)
##D plot(arbitrary.ocp)
##D points(occ,col="red",cex=0.5,pch=19)
##D 
##D ## with buffer:
##D mpa.ocp  <- ecospat.occupied.patch(pred.bin.mpa,occ, buffer=500000)
##D arbitrary.ocp  <- ecospat.occupied.patch(pred.bin.arbitrary,occ, buffer=500000)
##D 
##D plot(mpa.ocp) ## occupied patches: green area
##D points(occ,col="red",cex=0.5,pch=19)
##D plot(arbitrary.ocp)
##D points(occ,col="red",cex=0.5,pch=19)
## End(Not run)




cleanEx()
nameEx("ecospat.permut.glm")
### * ecospat.permut.glm

flush(stderr()); flush(stdout())

### Name: ecospat.permut.glm
### Title: GLM Permutation Function
### Aliases: ecospat.permut.glm

### ** Examples


## Not run: 
##D ecospat.permut.glm (get ("glm.Achillea_atrata", envir=ecospat.env), 1000)
## End(Not run)


cleanEx()
nameEx("ecospat.plot.kappa")
### * ecospat.plot.kappa

flush(stderr()); flush(stdout())

### Name: ecospat.plot.kappa
### Title: Plot Kappa
### Aliases: ecospat.plot.kappa
### Keywords: file

### ** Examples



Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
ecospat.plot.kappa(Pred, Sp.occ)



cleanEx()
nameEx("ecospat.plot.mess")
### * ecospat.plot.mess

flush(stderr()); flush(stdout())

### Name: ecospat.plot.mess
### Title: Plot MESS
### Aliases: ecospat.plot.mess

### ** Examples

## Not run: 
##D x <- ecospat.testData[c(2,3,4:8)]
##D proj <- x[1:90,] #A projection dataset.
##D cal <- x[91:300,] #A calibration dataset
##D 
##D #Create a MESS object 
##D mess.object <- ecospat.mess (proj, cal, w="default")
##D 
##D #Plot MESS 
##D ecospat.plot.mess (mess.object, cex=1, pch=15)
## End(Not run)



cleanEx()
nameEx("ecospat.plot.tss")
### * ecospat.plot.tss

flush(stderr()); flush(stdout())

### Name: ecospat.plot.tss
### Title: Plot True skill statistic (TSS)
### Aliases: ecospat.plot.tss
### Keywords: file

### ** Examples

Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
ecospat.plot.tss(Pred, Sp.occ)



cleanEx()
nameEx("ecospat.rand.pseudoabsences")
### * ecospat.rand.pseudoabsences

flush(stderr()); flush(stdout())

### Name: ecospat.rand.pseudoabsences
### Title: Sample Pseudo-Absences
### Aliases: ecospat.rand.pseudoabsences

### ** Examples

glob <- ecospat.testData[2:8]
presence <- ecospat.testData[c(2:3,9)]
presence <- presence[presence[,3]==1,1:2]
ecospat.rand.pseudoabsences (nbabsences=10, glob=glob, colxyglob=1:2, colvar = "all", 
presence= presence, colxypresence=1:2, mindist=20)



cleanEx()
nameEx("ecospat.rangesize")
### * ecospat.rangesize

flush(stderr()); flush(stdout())

### Name: ecospat.rangesize
### Title: Quantification of the range size of a species using habitat
###   suitability maps and IUCN criteria)
### Aliases: ecospat.rangesize
### Keywords: file

### ** Examples

 library(dismo)

# get predictor variables
fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
                     pattern='grd', full.names=TRUE )
predictors <- stack(fnames)
#plot(predictors)

# file with presence points
occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
colnames(occ) <- c("x","y")
occ <- ecospat.occ.desaggregation(occ,min.dist=1)

# fit a domain model, biome is a categorical variable
do <- domain(predictors, occ, factors='biome')

# predict to entire dataset
pred <- predict(do, predictors) 

plot(pred)
points(occ)


# use MPA to convert suitability to binary map
mpa.cutoff <- ecospat.mpa(pred,occ)

# use Boyce index to convert suitability to binary map
boyce <- ecospat.boyce(pred,  occ)
### use the boyce index to find a threshold
pred.bin.arbitrary <- ecospat.binary.model(pred,0.5)


pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)
names(pred.bin.mpa) <- "me.mpa"
pred.bin.arbitrary <- ecospat.binary.model(pred,0.5)
names(pred.bin.arbitrary) <- "me.arbitrary"

rangesize <- ecospat.rangesize(stack(pred.bin.mpa,pred.bin.arbitrary),
                               xy=occ,
                               resol=c(1,1),
                               eoo.around.modelocp =TRUE,
                               AOO.circles = TRUE,
                               d=200000,
                               lonlat =TRUE)


## Range size quantification
rangesize$RangeSize

names(rangesize$RangeObjects)


par(mfrow=c(1,3))

plot(ecospat.binary.model(pred,0),legend=FALSE, main="IUCN criteria")
## IUCN criteria & derivates
# plot AOO
plot(rangesize$RangeObjects$AOO,add=TRUE, col="red",legend=FALSE)
# plot EOO
plot(rangesize$RangeObjects$EOO@polygons,add=TRUE, border="red", lwd=2)
# plot circles around occurrences
plot(rangesize$RangeObjects$AOO.circle@polygons,add=TRUE,border="blue")

for(i in 1:2){
## plot the occupied patches of the model
plot(rangesize$RangeObjects$models.ocp[[i]],col=c("grey","blue","darkgreen"),
main=names(rangesize$RangeObjects$models.ocp[[i]]),legend=FALSE)
points(occ,col="red",cex=0.5,pch=19)
## plot EOO around model
plot(rangesize$RangeObjects$eoo.around.model[[i]]@polygons,add=TRUE,border="blue",lwd=2)
## plot EOO around occupied patches
plot(rangesize$RangeObjects$eoo.around.mo.ocp[[i]]@polygons,add=TRUE,border="darkgreen",
lwd=2)
## plot the modeled area within EOO
#plot(rangesize$RangeObjects$model.within.eoo[[i]],col=c("grey","blue","darkgreen"),legend=FALSE)
#points(occ,col="red",cex=0.5,pch=19)
#plot(rangesize$RangeObjects$EOO@polygons,add=TRUE, border="red", lwd=2)
}

### Alpha-hulls are not included in the function yet because of Licence limitations.
### However, alpha-hulls can easily be included manually (see also the help file of 
### the alpha hull package):

#  require(alphahull)
#  alpha = 2 # alpha value of 2 recommended by IUCN
  
#  del<-delvor(occ)
#  dv<-del$mesh
#  mn <- mean(sqrt(abs(del$mesh[,3]-del$mesh[,5])^2+abs(del$mesh[,4]-del$mesh[,6])^2))*alpha
#  alpha.hull<-ahull(del,alpha=mn) 
  
#  #Size of alpha-hulls
#  areaahull(h)


# plot alphahulls
# plot(rangesize$RangeObjects$models.ocp[[i]],col=c("grey","blue","darkgreen"),
#  main=names(rangesize$RangeObjects$models.ocp[[i]]),legend=FALSE)
# plot(alpha.hull,add=TRUE,lwd=1)
 
 


graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ecospat.rcls.grd")
### * ecospat.rcls.grd

flush(stderr()); flush(stdout())

### Name: ecospat.rcls.grd
### Title: Reclassifying grids function
### Aliases: ecospat.rcls.grd

### ** Examples


   ## Not run: 
##D bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
##D bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))
##D B3.rcl<-ecospat.rcls.grd(bio3,9) 
##D B12.rcl<-ecospat.rcls.grd(bio12,9)
##D B3B12.comb <- B12.rcl+B3.rcl*10
##D 
##D # Plotting a histogram of the classes
##D hist(B3B12.comb,breaks=100,col=heat.colors(88)) 
##D # Plotting the new RasterLayer (9x9 classes)
##D plot(B3B12.comb,col=rev(rainbow(88)),main="Stratified map") 
##D 
## End(Not run)


cleanEx()
nameEx("ecospat.recstrat_prop")
### * ecospat.recstrat_prop

flush(stderr()); flush(stdout())

### Name: ecospat.recstrat_prop
### Title: Random Ecologically Stratified Sampling of propotional numbers
### Aliases: ecospat.recstrat_prop

### ** Examples

  ## Not run: 
##D     bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
##D     bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))
##D     B3.rcl<-ecospat.rcls.grd(bio3,9) 
##D     B12.rcl<-ecospat.rcls.grd(bio12,9)
##D     B3B12.comb <- B12.rcl+B3.rcl*10
##D     B3B12.prop_samples <- ecospat.recstrat_prop(B3B12.comb,100)
##D     plot(B3B12.comb)
##D     points(B3B12.prop_samples$x,B3B12.prop_samples$y,pch=16,cex=0.6,col=B3B12.prop_samples$class)
##D     
##D     
##D   
## End(Not run)
  



cleanEx()
nameEx("ecospat.recstrat_regl")
### * ecospat.recstrat_regl

flush(stderr()); flush(stdout())

### Name: ecospat.recstrat_regl
### Title: Random Ecologically Stratified Sampling of equal numbers
### Aliases: ecospat.recstrat_regl

### ** Examples

  
  ## Not run: 
##D     bio3<- raster(system.file("external/bioclim/current/bio3.grd",package="biomod2"))
##D     bio12<- raster(system.file("external/bioclim/current/bio12.grd",package="biomod2"))
##D     B3.rcl<-ecospat.rcls.grd(bio3,9) 
##D     B12.rcl<-ecospat.rcls.grd(bio12,9)
##D     B3B12.comb <- B12.rcl+B3.rcl*10
##D     B3B12.regl_samples <- ecospat.recstrat_prop(B3B12.comb,100)
##D     plot(B3B12.comb)
##D     points(B3B12.regl_samples$x,B3B12.regl_samples$y,pch=16,cex=0.6,col=B3B12.regl_samples$class)
##D   
## End(Not run)




cleanEx()
nameEx("ecospat.sample.envar")
### * ecospat.sample.envar

flush(stderr()); flush(stdout())

### Name: ecospat.sample.envar
### Title: Sample Environmental Variables
### Aliases: ecospat.sample.envar

### ** Examples

## Not run: 
##D spp <- ecospat.testNiche
##D sp1 <- spp[1:32,1:3]
##D occ.sp1 <- ecospat.occ.desaggregation(dfvar=sp1,colxy=2:3,colvar=NULL, min.dist=500,plot=TRUE)
##D clim <- ecospat.testData[2:8]
##D 
##D occ_sp1 <- na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2,colspkept=1:2,dfvar=clim,
##D colvarxy=1:2,colvar="all",resolution=25))
## End(Not run)



cleanEx()
nameEx("ecospat.testData")
### * ecospat.testData

flush(stderr()); flush(stdout())

### Name: ecospat.testData
### Title: Test Data For The Ecospat package
### Aliases: ecospat.testData

### ** Examples

data(ecospat.testData)
str(ecospat.testData)
dim(ecospat.testData)
names(ecospat.testData)



cleanEx()
nameEx("ecospat.testEnvRaster")
### * ecospat.testEnvRaster

flush(stderr()); flush(stdout())

### Name: ecospat.testEnvRaster
### Title: Test Environmental Rasters for The Ecospat package
### Aliases: ecospat.testEnvRaster

### ** Examples

fpath <- system.file("extdata", "ecospat.testEnvRaster.RData", package="ecospat")
load(fpath)
plot(env)



cleanEx()
nameEx("ecospat.testMdr")
### * ecospat.testMdr

flush(stderr()); flush(stdout())

### Name: ecospat.testMdr
### Title: Test Data For The ecospat.mdr function
### Aliases: ecospat.testMdr

### ** Examples

data(ecospat.testMdr)
str(ecospat.testMdr)
dim(ecospat.testMdr)



cleanEx()
nameEx("ecospat.testNiche")
### * ecospat.testNiche

flush(stderr()); flush(stdout())

### Name: ecospat.testNiche
### Title: Test Data For The Niche Overlap Analysis
### Aliases: ecospat.testNiche

### ** Examples

data(ecospat.testNiche)
dim(ecospat.testNiche)
names(ecospat.testNiche)



cleanEx()
nameEx("ecospat.testNiche.inv")
### * ecospat.testNiche.inv

flush(stderr()); flush(stdout())

### Name: ecospat.testNiche.inv
### Title: Test Data For The Niche Dynamics Analysis In The Invaded Range
###   Of A Hypothetical Species
### Aliases: ecospat.testNiche.inv

### ** Examples

data(ecospat.testNiche.inv)
str(ecospat.testNiche.inv)
dim(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)



cleanEx()
nameEx("ecospat.testNiche.nat")
### * ecospat.testNiche.nat

flush(stderr()); flush(stdout())

### Name: ecospat.testNiche.nat
### Title: Test Data For The Niche Dynamics Analysis In The Native Range Of
###   A Hypothetical Species
### Aliases: ecospat.testNiche.nat

### ** Examples

data(ecospat.testNiche.nat)
str(ecospat.testNiche.nat)
dim(ecospat.testNiche.nat)
names(ecospat.testNiche.nat)



cleanEx()
nameEx("ecospat.testTree")
### * ecospat.testTree

flush(stderr()); flush(stdout())

### Name: ecospat.testTree
### Title: Test Tree For The Ecospat package
### Aliases: ecospat.testTree

### ** Examples

fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
tree <- read.tree(fpath)
plot(tree)



cleanEx()
nameEx("ecospat.varpart")
### * ecospat.varpart

flush(stderr()); flush(stdout())

### Name: ecospat.varpart
### Title: Variation Partitioning For GLM Or GAM
### Aliases: ecospat.varpart

### ** Examples

## Not run: 
##D ecospat.cv.example()
##D ecospat.varpart (model.1= get ("glm.Achillea_atrata", envir=ecospat.env), 
##D model.2= get ("glm.Achillea_millefolium", envir=ecospat.env), 
##D model.12= get ("glm.Achillea_millefolium", envir=ecospat.env))
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
