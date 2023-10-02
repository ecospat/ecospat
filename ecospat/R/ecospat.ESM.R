#####################################################################################
# Commented functions for calculating  "Ensemble of Small Models" (ESM) in R using
# the Biomod2 package
#
#
# These functions could be used to build ESMs using simple bivariate models which
# are averaged weighted by an evaluation score (e.g. AUC) accoring to
# Breiner et al (2015). They provide full functionality of the approach desribed in
# Breiner et al.
#
# a) ecospat.ESM.Modeling:
#    This function calibrates simple bivariate models as
#    in Lomba et al. 2010 and Breiner et al. 2015.
# b) ecospat.ESM.Projection:
#    Projects simple bivariate models on new.env
# c) ecospat.ESM.EnsembleModeling:
#    Evaluates and averages simple bivariate models by weighted means
#    to Ensemble Small Models as in Lomba et al. 2010 and Breiner et al. 2015.
# d) ecospat.ESM.EnsembleProjection:
#    Projecting calibrated ESMs into new space or time.
# e) ecospat.ESM.MergeModels:
#    Enables to merge different ecospat.ESM.Modeling outputs which were produced
#    e.g. on different computers or on a cluster
#
# Breiner et al. found for rare species ESMs to be on average superior to Maxent
# when compared by AUC or BOYCE. While BOYCE generally showed better performance of
# ESMs, AUC was in a few cases indifferent or decreased slightly
#
#

#
# REFERENCES:
# Breiner F.T., Guisan A., Bergamini A., Nobis M.P. (2015):
#  Overcoming limitations of modeling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218. doi: 10.1111/2041-210X.12403
# Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J. & Guisan, A. (2010):
#  Overcoming the rare species modeling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143, 2647-2657.
#
# Frank Breiner, Swiss Federal Research Institute WSL, frank.breiner@wsl.ch, Dec 2015
#####################################################################################



############################
# a)  ecospat.ESM.Modeling  		calibrates simple bivariate models
############################

#Description
# Calibrates simple bivariate models as in Lomba et al. 2010 and Breiner et al. 2015.

## FUNCTION'S ARGUMENTS
## data:                BIOMOD.formated.data object returned by BIOMOD_FormatingData
## NbRunEval:           number of dataset splitting replicates for the model evaluation (same as in Biomod)
## DataSplit:           percentage of dataset observations retained for the model training (same as in Biomod)
## DataSplitTable:      a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column corresponds to a 'RUN'. If filled, arguments NbRunEval, DataSplit and #do.full.models will be ignored.
## Prevalence:          either NULL or a 0-1 numeric used to build 'weighted response weights'. In contrast to Biomod the default is 0.5 (weighting presences equally to the absences). If NULL each observation (presence or absence) has the same weight (independent of the number of presences and absences).
## models:              vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF','MAXENT', "MAXNET" (same as in Biomod)
## modeling.id:	        character, the ID (=name) of modeling procedure. A random number by default.
## models.options:      BIOMOD.models.options object returned by BIOMOD_ModelingOptions (same as in Biomod). If none is provided standard ESM tuning parameterswill be used.
## tune:                logical. if true model tuning will be used to estimate optimal parameters for the models (Default: False).
## which.biva:          integer. which bivariate combinations should be used for modeling? Default: all
## weighting.score:     evaluation score used to weight single models to build ensembles: 'AUC', 'SomersD' (2xAUC-1), 'Kappa', 'TSS' or 'Boyce'
## ESM_Projection:      logical. set to FALSE, if no projections should be calculated. Only Evaluation scores based on test and train data will be returned. (Default: TRUE)
## new.env:             A set of explanatory variables onto which models will be projected . It could be a data.frame, a matrix, or a SpatRaster object. Make sure the column names (data.frame or matrix) or layer Names (SpatRaster) perfectly match with the names of variables used to build the models in the previous steps.
## parallel:            logical. If TRUE, the parallel computing is enabled (highly recommended)
## cleanup:             numeric. Not available. No cleanup is used by default.
## Yweights:            response points weights. This argument will only affect models that allow case weights. 

## Details:
# which biva allows to split model runs, e.g. if which.biva is 1:3, only the three first bivariate variable combinations will be modeled. This allows to run different biva splits on different computers.
# However, it is better not to use this option if all models should be run on a single computer.
# See ecospat_ESM_MergeModels for merging the single biva for ecospat_ESM_EnsembleModeling function. Default: running all biva models. NOTE: Make shure to give each of your biva runs a unique modeling.id!

## Values:
# modeling.id:    "character", id of modeling procedure
# models.:        location of models output
# models:         models argument
# pred.biva:      location of the bivariate models predictions
# calib.lines:    "BIOMOD.stored.array", calibrations and evaluation lines selected for each run
# NbRunEval:      Number of Evaluation run
# new.env:        set of explanatory variables onto which models were projected
# data:           'BIOMOD.formated.data' object
# wd:             working directory
# which.biva:     character. Which bivariate combinations should be modeled. If NULL (default) all possible combinations are considered.

## Authors:
# Frank Breiner <frank.breiner@unil.ch>
# with contributions from Mirko Di Febbraro
# with the updates of Flavien Collart

## References
#Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modeling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
# Breiner F.T., Guisan A., Bergamini A., Nobis M.P. (2015):
#  Overcoming limitations of modeling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218. doi: 10.1111/2041-210X.12403

##See Also
#ecospat.ESM.EnsembleModeling; ecospat.ESM.MergeModels

ecospat.ESM.Modeling <- function(data, NbRunEval = NULL, DataSplit = NULL, DataSplitTable = NULL, Prevalence = 0.5, weighting.score,
                                 models, tune = FALSE, modeling.id = as.character(format(Sys.time(), "%s")), models.options = NULL, which.biva = NULL,
                                 parallel = FALSE, cleanup = FALSE, Yweights = NULL) {
  if (!weighting.score %in% c("AUC", "TSS", "Boyce", "Kappa", 
                              "SomersD")) {
    stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")
  }
  if (data@has.data.eval) {
    stop("Evaluation with independant data is not supported yet!")
  }
  if ("PA.table" %in% slotNames(data)) {
    if (ncol(data@PA.table) > 1) {
      stop("It is not possible to use more than one Pseudo Absences dataset")
    }
    if (sum(data@PA.table[, 1]) != length(data@data.species)) {
      stop("The number of Pseudo Absences given in the argument resp.var of BIOMOD_FormatingData should be exactly the same as the in PA.nb.absences")
    }
  }
  if (weighting.score == "AUC" | weighting.score == "SomersD") {
    models.eval.meth <- "ROC"
  }
  if (weighting.score == "Kappa") {
    models.eval.meth <- "KAPPA"
  }
  if (weighting.score == "TSS") {
    models.eval.meth <- "TSS"
  }
  if (weighting.score == "Boyce") {
    models.eval.meth <- "ROC"
  }
  if (is.null(models.options)) {
    models.options <- biomod2::BIOMOD_ModelingOptions()
    models.options@GBM$n.trees <- 1000
    models.options@GBM$interaction.depth <- 4
    models.options@GBM$shrinkage <- 0.005
    models.options@GAM$select <- TRUE
    models.options@CTA$control$cp <- 0
    models.options@ANN$size <- 8
    models.options@ANN$decay <- 0.001
    models.options@MARS$interaction.level <- 0
    models.options@MARS$nprune <- 2
    models.options@MAXENT$product <- FALSE
    models.options@MAXENT$threshold <- FALSE
    models.options@MAXENT$betamultiplier <- 0.5
    models.options@GLM$test <- "none"
  }
  if ("MAXENT" %in% models) {
    if (!file.exists(paste(models.options@MAXENT$path_to_maxent.jar, 
                           "maxent.jar", sep = "/"))) {
      stop("maxent.jar file not found!")
    }
  }
  if (is.null(NbRunEval) & is.null(DataSplitTable)) {
    stop("Need to give a value for NbRunEval  yhen DataSplitTable is null")
  }
  if (is.null(DataSplit) & is.null(DataSplitTable)) {
    stop("Need to give a value for DataSplit  when DataSplitTable is null")
  }
  models <- sort(models)
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  dir.create(paste("./ESM.BIOMOD.output", data@sp.name, sep = "_"))
  newwd <- paste(getwd(), "/ESM.BIOMOD.output_", data@sp.name, 
                 sep = "")
  setwd(newwd)
  combinations <- combn(colnames(data@data.env.var), 2)
  if (is.null(which.biva)) {
    which.biva <- 1:ncol(combinations)
  }
  mydata <- data
  if (is.null(DataSplitTable)) {
    calib.lines <- .CreatingDataSplitTable(bm.format = mydata, 
                                                     NbRunEval = NbRunEval, DataSplit = DataSplit)
    if ("PA.table" %in% slotNames(data)){
      
      colnames(calib.lines) = paste0("_PA1",colnames(calib.lines))
      
    }else{
      colnames(calib.lines) = paste0("_allData",colnames(calib.lines))
      
    }
    
    
  }
  else {
    calib.lines <- DataSplitTable
    calib.lines <- cbind(calib.lines,TRUE)
    if ("PA.table" %in% slotNames(data)){
      
      colnames(calib.lines)[ncol(calib.lines)]<- "_PA1_Full"
      
    }else{
      colnames(calib.lines)[ncol(calib.lines)]<-"_allData_Full"
      
    }
  }
  if (is.null(NbRunEval)) {
    if (ncol(calib.lines) > 1) {
      if (sum(!as.data.frame(calib.lines)[, ncol(calib.lines)]) == 
          0) {
        NbRunEval <- ncol(calib.lines) - 1
      }
      else {
        NbRunEval <- ncol(calib.lines)
      }
    }
    else {
      if (sum(!calib.lines) == 0) {
        NbRunEval <- ncol(calib.lines) - 1
      }
      else {
        NbRunEval <- ncol(calib.lines)
      }
    }
  }
  mymodels <- list()
  if (parallel == FALSE) {
    for (k in which.biva) {
      mydata@data.env.var <- data@data.env.var[, colnames(data@data.env.var) %in% 
                                                 combinations[, k]]
      mydata@sp.name <- paste("ESM.BIOMOD", k, sep = ".")
      if (tune == TRUE) {
        models.options <- biomod2::BIOMOD_Tuning(bm.format = mydata, 
                                                 models = models[models != "RF"], bm.options = models.options, 
                                                 weights = Yweights)$models.options
      }
      mymodels[[k]] <- "failed"
      try(mymodels[[k]] <- biomod2::BIOMOD_Modeling(bm.format = mydata, 
                                                    models = models, bm.options = models.options, 
                                                    CV.strategy="user.defined",
                                                    CV.nb.rep = NbRunEval, metric.eval = models.eval.meth, 
                                                    CV.user.table = as.matrix(calib.lines[,-ncol(calib.lines)]), prevalence = Prevalence, 
                                                    CV.do.full.models = TRUE, var.import = 0, modeling.id = modeling.id, 
                                                    weights = Yweights))
      if (cleanup != FALSE) {
        warning("Unfortunately cleanup is no more available. Please use terra::tmpFiles after this function to remove temporary files.")
      }
    }
  }
  if (parallel == TRUE) {
    mymodels <- foreach::foreach(k = which.biva, .packages = c("biomod2", 
                                                      "terra")) %dopar% {
                                                        setwd(newwd)
                                                        mydata@data.env.var <- data@data.env.var[, colnames(data@data.env.var) %in% 
                                                                                                   combinations[, k]]
                                                        mydata@sp.name <- paste("ESM.BIOMOD", k, sep = ".")
                                                        if (cleanup != FALSE) {
                                                          warning("Unfortunately cleanup is no more available. Please use terra::tmpFiles after this function to remove temporary files.")

                                                        }
                                                        if (tune == TRUE) {
                                                          models.options <- biomod2::BIOMOD_Tuning(bm.format = mydata, 
                                                                                                   models = models[models != "RF"], bm.options = models.options, 
                                                                                                   weights = Yweights)$models.options
                                                        }
                                                        biomod2::BIOMOD_Modeling(bm.format = mydata, 
                                                                                 models = models, bm.options = models.options, 
                                                                                 CV.strategy="user.defined",
                                                                                 CV.nb.rep = NbRunEval, metric.eval = models.eval.meth, 
                                                                                 CV.user.table = as.matrix(calib.lines[,-ncol(calib.lines)]), prevalence = Prevalence, 
                                                                                 CV.do.full.models = TRUE, var.import = 0, modeling.id = modeling.id, 
                                                                                 weights = Yweights)
                                                      }
  }
  output <- list(modeling.id = modeling.id, models. = grep(modeling.id, 
                                                           gtools::mixedsort(list.files(getwd(), "models.out", recursive = TRUE, 
                                                                                        full.names = TRUE)), value = TRUE), models = models, 
                 calib.lines = calib.lines, NbRunEval = NbRunEval, data = data, 
                 wd = getwd(), which.biva = which.biva, mymodels = mymodels)
  save(output, file = paste("ESM_Modeling..models", modeling.id, 
                            "out", sep = "."))
  setwd(iniwd)
  return(output)
}

############################
# b)  ecospat.ESM.Projection    projects simple bivariate models on new.env.
############################

## FUNCTION'S ARGUMENTS
## ESM.modeling.output:   BIOMOD.formated.data object returned by BIOMOD_FormatingData
## new.env:               A set of explanatory variables onto which models will be projected . It could be a data.frame, a matrix, and a SpatRaster object. Make sure the column names (data.frame or matrix) or layer Names (SpatRaster) perfectly match with the names of variables used to build the models in the previous steps.
## parallel:              logical. If TRUE, the parallel computing is enabled
## cleanup:               numeric. Not available. No cleanup is used by default.

## Details:
# The basic idea of ensemble of small models (ESMs) is to model a species distribution based on
# small, simple models, for example all possible bivariate models (i.e.  models that contain only two
# predictors at a time out of a larger set of predictors), and then combine all possible bivariate models
# into an ensemble (Lomba et al. 2010; Breiner et al. 2015).
# The ESM set of functions could be used to build ESMs using simple bivariate models which are
# averaged using weights based on model performances (e.g. AUC) accoring to Breiner et al (2015).
# They provide full functionality of the approach described in Breiner et al. (2015).
# The name of new.env must be a regular expression (see ?regex)

## Values:
# Returns the projections for all selected models (same as in biomod2) See "BIOMOD.projection.out" for details.

## Authors:
# Frank Breiner <frank.breiner@unil.ch>
# with the updates of Flavien Collart

## References
# Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modeling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
# Breiner F.T., Guisan A., Bergamini A., Nobis M.P. (2015): Overcoming limitations of modeling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218. doi: 10.1111/2041-210X.12403

##See Also
# ecospat.ESM.Modeling; ecospat.ESM.EnsembleModeling; ecospat.ESM.EnsembleProjection


ecospat.ESM.Projection <- function(ESM.modeling.output, new.env, name.env = NULL, parallel = FALSE, cleanup = FALSE) {
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(ESM.modeling.output$wd)
  models <- ESM.modeling.output$models
  models. <- ESM.modeling.output$models.
  mymodels <- ESM.modeling.output$mymodels
  data <- ESM.modeling.output$data
  combinations <- combn(colnames(ESM.modeling.output$data@data.env.var), 
                        2)
  which.biva <- ESM.modeling.output$which.biva
  NbRunEval <- ESM.modeling.output$NbRunEval
  modeling.id <- ESM.modeling.output$modeling.id
  if (is.matrix(new.env)) {
    new.env <- as.data.frame(new.env)
  }
  if (is.null(name.env)) 
    name.env <- deparse(substitute(new.env))
  if (parallel == FALSE) {
    
    for (k in 1:length(mymodels)) {
      mymodel <- mymodels[[k]]
      dir.create(path = paste0("./","ESM.BIOMOD.",k,"/proj_",paste(name.env, "ESM.BIOMOD", k, modeling.id, sep = ".")))
      
      if (is.character(mymodel)) {
        (next)()
      }
      if (is.data.frame(new.env)) {
        
        newdata=new.env[, colnames(new.env) %in% combinations[, k]]
        if("PA.table" %in% slotNames(data)){
          
        models.chosen = grep("allRun", 
                             mymodel@models.computed[-grep("allData",mymodel@models.computed)], 
                             value = TRUE)}
        else{
          models.chosen = grep("allRun",
                               mymodel@models.computed, 
                                                    value = TRUE)
        }
                             
        for(i in 1:length(models)){ ## Error with BIOMOD_Projection function
          modelToProject <- get(biomod2::BIOMOD_LoadModels(mymodel,full.name = grep(paste0("_",models[i]),models.chosen,value = TRUE)))
          if(models[i]=="MAXENT"){
            map <- biomod2::predict(modelToProject,newdata=newdata,temp_workdir = modelToProject@model_output_dir,overwrite=T,on_0_1000 = TRUE, omit.na = TRUE)
          }else{
            map <- biomod2::predict(modelToProject,newdata=newdata,overwrite=T,on_0_1000 = TRUE, omit.na = TRUE)
          }
          
          save(map,file=file.path(paste0("./","ESM.BIOMOD.",k,"/proj_",paste(name.env, "ESM.BIOMOD", k, modeling.id, sep = "."),"/proj_",paste(name.env,models[i] ,"ESM.BIOMOD", k, modeling.id, sep = "."),".RData")))
          
        
        }
      }
      if (inherits(new.env, "SpatRaster")) {
        
        newdata=newdata=subset(new.env,combinations[, k])
        if("PA.table" %in% slotNames(data)){
          
          models.chosen = grep("allRun", 
                               mymodel@models.computed[-grep("allData",mymodel@models.computed)], 
                               value = TRUE)}
        else{
          models.chosen = grep("allRun",
                               mymodel@models.computed, 
                               value = TRUE)
        }
        
        for(i in 1:length(models)){
          modelToProject <- get(biomod2::BIOMOD_LoadModels(mymodel,full.name = grep(paste0("_",models[i]),models.chosen,value = TRUE)))
          if(models[i]=="MAXENT"){
            map <- biomod2::predict(modelToProject,newdata=newdata,temp_workdir = modelToProject@model_output_dir,overwrite=T,on_0_1000 = TRUE, omit.na = TRUE)
          }else{
            map <- biomod2::predict(modelToProject,newdata=newdata,overwrite=T,on_0_1000 = TRUE, omit.na = TRUE)
          }
          
          terra::writeRaster(map,paste0("./","ESM.BIOMOD.",k,"/proj_",paste(name.env, "ESM.BIOMOD", k, modeling.id, sep = "."),"/proj_",paste(name.env,models[i] ,"ESM.BIOMOD", k, modeling.id, sep = "."),".tif"),overwrite =TRUE)
        }
        
        
      }      
      }
    }
  if (parallel == TRUE) {
    
    if (inherits(new.env, "SpatRaster")) {
      new.env <- terra::wrap(new.env) ## Allow parallelisation
    }
    
    for(g in 1:length(mymodels)){
      dir.create(path = paste0("./","ESM.BIOMOD.",g,"/proj_",paste(name.env, "ESM.BIOMOD", g, modeling.id, sep = ".")))
    }
    
    foreach::foreach(k = 1:length(mymodels), .packages = c("biomod2",
                                                  "terra","base")) %dopar% {
                                                    mymodel <- mymodels[[k]]
                                                    
                                                    if (!(is.character(mymodel))) {
                                                      i=0
                                                      setwd(ESM.modeling.output$wd)

                                                      if (is.data.frame(new.env)) {
                                                        
                                                        newdata=new.env[, colnames(new.env) %in% combinations[, k]]
                                                        if("PA.table" %in% slotNames(data)){
                                                          
                                                          models.chosen = grep("allRun", 
                                                                               mymodel@models.computed[-grep("allData",mymodel@models.computed)], 
                                                                               value = TRUE)}
                                                        else{
                                                          models.chosen = grep("allRun",
                                                                               mymodel@models.computed, 
                                                                               value = TRUE)
                                                        }
                                                        
                                                        for(i in 1:length(models)){
                                                          modelToProject <- get(biomod2::BIOMOD_LoadModels(mymodel,full.name = grep(paste0("_",models[i]),models.chosen,value = TRUE)))
                                                          if(models[i]=="MAXENT"){
                                                            map <- biomod2::predict(modelToProject,newdata=newdata,temp_workdir = modelToProject@model_output_dir,overwrite=T,on_0_1000 = TRUE, omit.na = TRUE)
                                                          }else{
                                                            map <- biomod2::predict(modelToProject,newdata=newdata,overwrite=T,on_0_1000 = TRUE, omit.na = TRUE)
                                                          }
                                                          
                                                          save(map,file=file.path(paste0("./","ESM.BIOMOD.",k,"/proj_",paste(name.env, "ESM.BIOMOD", k, modeling.id, sep = "."),"/proj_",paste(name.env,models[i] ,"ESM.BIOMOD", k, modeling.id, sep = "."),".RData")))
                                                          
                                                          
                                                        }
                                                      }else{
                                                          new.env <- terra::rast(new.env)
                                                          
                                                          newdata <- terra::subset(new.env, combinations[, k])
                                                          new.env <- terra::wrap(new.env)
                                                          if("PA.table" %in% slotNames(data)){
                                                            
                                                            models.chosen = grep("allRun", 
                                                                                 mymodel@models.computed[-grep("allData",mymodel@models.computed)], 
                                                                                 value = TRUE)}
                                                          else{
                                                            models.chosen = grep("allRun",
                                                                                 mymodel@models.computed, 
                                                                                 value = TRUE)
                                                          }
                                                          
                                                          for(i in 1:length(models)){
                                                            modelToProject <- get(biomod2::BIOMOD_LoadModels(mymodel,full.name = grep(paste0("_",models[i]),models.chosen,value = TRUE)))
                                                            if(models[i]=="MAXENT"){
                                                              map <- biomod2::predict(modelToProject,newdata=newdata,temp_workdir = modelToProject@model_output_dir,on_0_1000 = TRUE, omit.na = TRUE,overwrite=T)
                                                            }else{
                                                              map <- biomod2::predict(modelToProject,newdata=newdata,on_0_1000 = TRUE, omit.na = TRUE,overwrite=T)
                                                            }
                                                            
                                                            terra::writeRaster(map,file=file.path(paste0(ESM.modeling.output$wd,"/","ESM.BIOMOD.",k,"/proj_",paste(name.env, "ESM.BIOMOD", k, modeling.id, sep = "."),"/proj_",paste(name.env,models[i] ,"ESM.BIOMOD", k, modeling.id, sep = "."),".tif")),overwrite=T)
                                                            
                                                            
                                                          }
                                                        
                                                        }
    
                                                    }
                                                  }
    if(!is.data.frame(new.env)){
      new.env <- terra::rast(new.env)
    }
                                                  
  }
  if (cleanup != FALSE) {
    warning("Unfortunately cleanup is no more available. Please use terra::tmpFiles after this function to remove temporary files.")
  }
  output <- list(proj.name = name.env, modeling.id = modeling.id, 
                 models. = grep(modeling.id, gtools::mixedsort(list.files(getwd(), 
                                                                          "models.out", recursive = TRUE, full.names = TRUE)), 
                                value = TRUE), models = models, pred.biva = grep(modeling.id, 
                                                                                 gtools::mixedsort(list.files(getwd(), paste("proj_", 
                                                                                                                             name.env, sep = ""), recursive = TRUE, full.names = TRUE)), 
                                                                                 value = TRUE), NbRunEval = NbRunEval, name.env = name.env, 
                 new.env.raster = inherits(new.env, "SpatRaster"), wd = getwd(), 
                 which.biva = which.biva)
  save(output, file = paste("ESM_Projections", name.env, modeling.id, 
                            "out", sep = "."))
  setwd(iniwd)
  return(output)
}



############################
# c)  ecospat.ESM.EnsembleModeling		Evaluates and averages simple bivariate models to	ESMs
############################


#Description
## Evaluates and averages simple bivariate models  by weighted means to Ensemble Small Models as in Lomba et al. 2010 and Breiner et al. 2015.

## FUNCTION'S ARGUMENTS
## ESM.modeling.output:   BIOMOD.formated.data object returned by BIOMOD_FormatingData
## weighting.score:       evaluation score used to weight single models to build ensembles: 'AUC', 'SomersD' (2xAUC-1), 'Kappa', 'TSS' or 'Boyce'
## threshold:             models evaluated equal or worse than the threshold will be removed from building the ensembles (Default: 0.5 AUC and 0 for the other inidces; a threshold of 0.5 for AUC and 0 for SomersD and Boyce index will remove models no better than random from the ensembles)

## Values:
# species:          species name
# ESM.fit:          data.frame of the predicted values for the data used to build the models.
# ESM.evaluations:  data.frame with evaluations scores for the ESMs
# weights:          weighting scores used to weight the bivariate models to build the single ESM
# weights.EF:       weighting scores used to weight the single ESM to build the ensemble of ESMs from different modelling techniques (only available if >1 modelling techniques were selected).
# failed:           bivariate models which failed because they could not be calibrated.

## Authors:
# Frank Breiner
# With the updates of Flavien Collart

## References
#Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
#Breiner, F. B., A. Guisan, A. Bergamini, M. P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, submitted.

## See also
# ecospat.ESM.Modeling; ecospat.ESM.MergeModels

ecospat.ESM.EnsembleModeling <- function(ESM.modeling.output, weighting.score, threshold = NULL, models){
  if (!weighting.score %in% c("AUC", "TSS", "Boyce", "Kappa", 
                              "SomersD")) {
    stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")
  }
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(ESM.modeling.output$wd)
  data <- ESM.modeling.output$data
  models. <- ESM.modeling.output$models.
  NbRunEval <- ESM.modeling.output$NbRunEval
  models <- ESM.modeling.output$models
  calib.lines <- as.data.frame(ESM.modeling.output$calib.lines)
  ESM_Projection <- ESM.modeling.output$ESM_Projection
  new.env <- ESM.modeling.output$new.env
  for (i in 1:length(models.)) load(models.[i])
  models. <- NULL
  for (n in 1:length(ESM.modeling.output$modeling.id)) {
    tmp.list.files <- grep(ESM.modeling.output$modeling.id[n], 
                           grep(paste(".", ESM.modeling.output$modeling.id[n], 
                                      ".models.out", sep = ""), ls(), value = TRUE, 
                                fixed = TRUE), value = TRUE)
    models. <- c(models., tmp.list.files[gtools::mixedorder(gsub(".", 
                                                         "_", tmp.list.files, fixed = TRUE))])
  }
  mymodel <- list()
  for (i in 1:length(models.)) mymodel[[i]] <- get(models.[i])
  weights <- unlist(lapply(mymodel, function(x) {
    y <- x
    if (weighting.score == "Boyce") {
    ## Possible errors can appear in this section
      z <- biomod2::get_predictions(y, model.as.col = TRUE)
      z <- z[, -grep("allRun",colnames(z))]
      y.eval <- biomod2::get_evaluations(y)
      if("PA.table" %in% slotNames(data)){
        y.eval <- y.eval[y.eval$PA != "allData",] 
        y.eval$run[y.eval$run=="allRun"] = "Full"
      }else{
        y.eval$run[y.eval$run=="allRun"] = "Full"
      }

      x <- matrix(y.eval$validation, nrow = length(unique(y.eval$algo)), 
                  ncol = NbRunEval + 1, byrow = FALSE)
      colnames(x) <- unique(y.eval$run)
      rownames(x) <- unique(y.eval$algo)
      
      if(anyNA(data@data.species)){
        data@data.species[is.na(data@data.species)] =0
      }
      
      if (length(models) > 1) {
        x[, ] <- NA
        x <- x[, colnames(x) != "Full" & colnames(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
      }
      else {
        x[] <- NA
        x <- x[, colnames(x) != "Full" & colnames(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
      }
      for (n in 1:(length(models))) {
        model <- models[n]
        if (!TRUE %in% c(grepl("Full", y@models.failed), 
                         grepl(paste("RUN", NbRunEval + 1, sep = ""), 
                               y@models.failed))) {
          for (i in 1:NbRunEval) {
            if (sum(is.na(z[, grep(paste("RUN", i, "_", 
                                         model, sep = ""), colnames(z))])) != nrow(z)) {
              if (length(models) > 1) {
                
                  x[rownames(x) == model, colnames(x) == 
                      paste("RUN", i, sep = "")] <- ecospat.boyce(z[!calib.lines[,i ],grep(paste("RUN", i, "_", model, sep = ""), colnames(z))], 
                                                                  z[ data@data.species ==  1 & !calib.lines[,i],grep(paste("RUN", i, "_", model, 
                                                                                                                           sep = ""), colnames(z))], PEplot = F)$cor
                

              }else {

                  x[names(x) == paste("RUN", i, sep = "")] <- ecospat.boyce(z[!calib.lines[,i ],grep(paste("RUN", i, "_", model, sep = ""), colnames(z))], 
                                                                            z[data@data.species ==  1 & ! calib.lines[,i],grep(paste("RUN", i, "_", model, 
                                                                                                                                     sep = ""), colnames(z))], PEplot = F)$cor
                
              }
            }
          }
        }
      }
      if (length(models) > 1) {
        x <- round(apply(x, 1, mean, na.rm = TRUE), 4)
      }
      else {
        x <- round(mean(x, na.rm = TRUE), 4)
      }
    }
    else {
      y.eval <- biomod2::get_evaluations(y)
      if("PA.table" %in% slotNames(data)){
        y.eval <- y.eval[y.eval$PA != "allData",] 
        y.eval$run[y.eval$run=="allRun"] = "Full"
      }else{
        y.eval$run[y.eval$run=="allRun"] = "Full"
      }
      
      x <- matrix(y.eval$validation, nrow = length(unique(y.eval$algo)), 
                  ncol = NbRunEval + 1, byrow = FALSE)
      x.calib <- matrix(y.eval$calibration, nrow = length(unique(y.eval$algo)), 
                        ncol = NbRunEval + 1, byrow = FALSE)
      colnames(x) <- unique(y.eval$run)
      rownames(x) <- unique(y.eval$algo)
      colnames(x.calib) <- unique(y.eval$run)
      rownames(x.calib) <- unique(y.eval$algo)
      if (length(models) > 1) {
        for (row in models) {
          if (ncol(calib.lines) == NbRunEval + 1) {
            if (is.na(x.calib[row, NbRunEval + 1])) {
              x[row, ] <- NA
            }
          }
        }
        x <- x[, colnames(x) != "Full" & colnames(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
        if (NbRunEval > 1) {
          x <- round(apply(x, 1, mean, na.rm = TRUE), 
                     4)
        }
      }
      else {
        if (ncol(calib.lines) == NbRunEval + 1) {
          if (is.na(x.calib[NbRunEval + 1])) {
            x <- NA
          }
        }
        x <- x[, colnames(x) != "Full" & colnames(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
        x <- round(mean(x, na.rm = TRUE), 4)
      }
      if (weighting.score == "SomersD") {
        x <- x * 2 - 1
      }
    }
    if(length(models)==1){
      names(x) <- paste(unique(y.eval$algo), y@sp.name, sep = ".")
    }else{ ## Here need to test
      names(x) <- paste(names(x), y@sp.name, sep = ".")
    }
    
    return(x)
  }), recursive = TRUE)
  failed.mod <- "none"
  failed.mod <- lapply(mymodel, function(x) {
    return(if (x@models.failed[1] == "none") {
      paste(x@sp.name, "none", sep = ": ")
    } else {
      x@models.failed
    })
  })
  failed <- "none"
  if (sum(is.na(weights)) > 0) {
    warning(cat(paste("\n\n\n################################\n", 
                      paste("The following bivariate model(s) failed and are not included for building the ESMs\n", 
                            sep = " "))), print(names(which(is.na(weights)))))
    failed <- names(which(is.na(weights)))
  }
  if (is.null(threshold)) {
    if (weighting.score == "AUC") {
      threshold <- 0.5
    }
    else {
      threshold <- 0
    }
  }
  if (sum(weights <= threshold & !is.na(weights)) > 0) {
    warning(cat(paste("\n\n\n################################\n", 
                      paste("The following bivariate model(s) is (are) not included for building the ESMs\n because evaluation score is smaller or equal to the given threshold of", 
                            weighting.score, "=", threshold, ":\n", sep = " "))), 
            print(names(weights[weights <= threshold & !is.na(weights)])))
  }
  weights[(weights <= threshold) | is.na(weights)] <- 0
  if (length(models) == 1) {
    mymodel <- mymodel[which(weights > 0)]
    weights <- weights[which(weights > 0)]
  }
  test.pred <- lapply(mymodel, function(x) {
    x <- biomod2::get_predictions(x, model.as.col = TRUE)
    return(x)
  })
  test.ESM <- NULL
  biva.st2 <- do.call(cbind, test.pred)
  for (i in 1:length(models)) {
    for (run in 1:(ncol(calib.lines)-1)) {
      Model.biva.prediction <- biva.st2[, grep(paste("RUN", 
                                                     run, "_", models[i], sep = ""), colnames(biva.st2))]
      ModelToWeight <- paste(models[i], sub("_.*", "", 
                                            colnames(Model.biva.prediction)), sep = ".")
      weights2 <- weights[ModelToWeight]
      test.ESM1 <- apply(Model.biva.prediction, 1, function(x) weighted.mean(x, 
                                                                             weights2, na.rm = TRUE))
      test.ESM <- cbind(test.ESM, test.ESM1)
      colnames(test.ESM)[ncol(test.ESM)] <- paste("RUN", 
                                                  run, "_", models[i], sep = "")
    }
    
    ## FullModel
    Model.biva.prediction <- biva.st2[, grep(paste("allRun_", models[i], sep = ""), colnames(biva.st2))]
    ModelToWeight <- paste(models[i], sub("_.*", "", 
                                          colnames(Model.biva.prediction)), sep = ".")
    weights2 <- weights[ModelToWeight]
    test.ESM1 <- apply(Model.biva.prediction, 1, function(x) weighted.mean(x, 
                                                                           weights2, na.rm = TRUE))
    test.ESM <- cbind(test.ESM, test.ESM1)
    colnames(test.ESM)[ncol(test.ESM)] <- paste("Full_", models[i], sep = "")
  }
  test.ESM <- as.data.frame(test.ESM)
  DATA <- cbind(1:length(data@data.species), resp.var = data@data.species, 
                test.ESM/1000)
  if ("PA.table" %in% slotNames(data)) {
    DATA$resp.var[is.na(DATA$resp.var)] <- 0
  }
  EVAL <- NULL
  for (i in 1:NbRunEval) {
    DATA1 <- cbind(DATA[, 1:2], DATA[, grep(paste("RUN", 
                                                  i, "_", sep = ""), colnames(DATA))])
    colnames(DATA1)[3] <- colnames(DATA)[2 + i]
    if (length(models) > 1) {
      if (nrow(DATA1) %in% apply(DATA1[, 3:ncol(DATA1)], 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        failed <- which(nrow(DATA1) == apply(DATA1, 2, 
                                             function(x) {
                                               sum(is.na(x))
                                             }))
        DATA1[, nrow(DATA1) == apply(DATA1, 2, function(x) {
          sum(is.na(x))
        })] <- 0
      }
    }
    else {
      if (nrow(DATA1) == sum(is.na(DATA1[, 3]))) {
        failed <- paste("RUN", i, sep = "")
        DATA1[, 3] <- 0
      }
    }
    EVAL1 <- PresenceAbsence::presence.absence.accuracy(DATA1[!calib.lines[, 
                                                                           i], ], threshold = as.vector(PresenceAbsence::optimal.thresholds(DATA1[!calib.lines[, 
                                                                                                                                                               i], ], opt.methods = "MaxSens+Spec")[-1], mode = "numeric"))
    EVAL1 <- EVAL1[c(1, 2, 4:7, 9:12)]
    EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 
      1
    if (length(models) > 1) {
      if (nrow(DATA1) %in% apply(DATA1[, 3:ncol(DATA1)], 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        EVAL1[EVAL1$model == names(failed), 2:11] <- NA
      }
    }
    else {
      if (nrow(DATA1) == sum(is.na(DATA1[, 3]))) {
        EVAL1[EVAL1$model == names(failed), 2:11] <- NA
      }
    }
    EVAL1$SomersD <- EVAL1$AUC * 2 - 1
    EVAL1$Boyce <- EVAL1$MPA <- NA
    for (n in 1:nrow(EVAL1)) {
      EVAL1$MPA[EVAL1$model == EVAL1$model[n]] <- ecospat.mpa(DATA1[!calib.lines[, 
                                                                                 i] & DATA1[, 2] == 1, EVAL1$model[n]])
      EVAL1$Boyce[EVAL1$model == EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[, 
                                                                                     i], EVAL1$model[n]], DATA1[!calib.lines[, i] & 
                                                                                                                  DATA1[, 2] == 1, EVAL1$model[n]], PEplot = F)$cor
    }
    EVAL1$technique <- unlist(strsplit(EVAL1$model, split = "_"))[seq(2, 
                                                                      nrow(EVAL1) * 2, 2)]
    EVAL1$RUN <- paste("RUN", i, sep = "")
    EVAL <- rbind(EVAL, EVAL1)
  }
  if (length(models) > 1) {
    weights.double <- stats::aggregate(EVAL[, weighting.score], 
                                by = list(EVAL$technique), FUN = mean, na.rm = TRUE)
    weights.double <- weights.double[order(weights.double[, 
                                                          1]), ]
    if (!is.null(threshold)) {
      weights.double[weights.double <= threshold] <- 0
      if (sum(weights.double == 0) > 0) {
        warning(cat(paste("\n\n\n################################\n", 
                          paste("The following ESM model(s) is (are) not included for building the 'double ensemble'\n because evaluation score is smaller or equal to the given threshold:\n", 
                                list(weights.double[weights.double[, 2] == 
                                                      0, 1]), sep = " ", "\n################################\n"))))
      }
    }
    for (run in 1:(ncol(calib.lines)-1)) {
      test.ESM$EF <- apply(test.ESM[, grep(paste("RUN", 
                                                 run, "_", sep = ""), colnames(test.ESM))], 1, 
                           function(x) weighted.mean(x, weights.double[, 
                                                                       2], na.rm = TRUE))
      colnames(test.ESM)[ncol(test.ESM)] <- paste("RUN", 
                                                  run, "_EF", sep = "")
    }
    DATA <- cbind(1:length(data@data.species), resp.var = data@data.species, 
                  test.ESM/1000)
    if ("PA.table" %in% slotNames(data)) {
      DATA$resp.var[is.na(DATA$resp.var)] <- 0
    }
    EVAL <- NULL
    for (i in 1:NbRunEval) {
      DATA1 <- cbind(DATA[, 1:2], DATA[, grep(paste("RUN", 
                                                    i, "_", sep = ""), colnames(DATA))])
      if (nrow(DATA1) %in% apply(as.data.frame(DATA1[, 3:ncol(DATA1)]), 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        failed <- which(nrow(DATA1) == apply(DATA1, 2, 
                                             function(x) {
                                               sum(is.na(x))
                                             }))
        DATA1[, nrow(DATA1) == apply(DATA1, 2, function(x) {
          sum(is.na(x))
        })] <- 0
      }
      EVAL1 <- PresenceAbsence::presence.absence.accuracy(DATA1[!calib.lines[, 
                                                            i], ], threshold = as.vector(PresenceAbsence::optimal.thresholds(DATA1[!calib.lines[, 
                                                                                                                               i], ], opt.methods = "MaxSens+Spec")[-1], mode = "numeric"))
      EVAL1 <- EVAL1[c(1, 2, 4:7, 9:12)]
      EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 
        1
      if (nrow(DATA1) %in% apply(as.data.frame(DATA1[, 3:ncol(DATA1)]), 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        EVAL1[EVAL1$model == names(failed), 2:11] <- NA
      }
      EVAL1$SomersD <- EVAL1$AUC * 2 - 1
      EVAL1$Boyce <- EVAL1$MPA <- NA
      for (n in 1:nrow(EVAL1)) {
        EVAL1$MPA[EVAL1$model == EVAL1$model[n]] <- ecospat.mpa(DATA1[!calib.lines[, 
                                                                                   i] & DATA1[, 2] == 1, EVAL1$model[n]])
        EVAL1$Boyce[EVAL1$model == EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[, 
                                                                                       i], EVAL1$model[n]], DATA1[!calib.lines[, i] & 
                                                                                                                    DATA1[, 2] == 1, EVAL1$model[n]], PEplot = F)$cor
      }
      EVAL1$technique <- unlist(strsplit(EVAL1$model, split = "_"))[seq(2, 
                                                                        nrow(EVAL1) * 2, 2)]
      EVAL1$RUN <- paste("RUN", i, sep = "")
      EVAL <- rbind(EVAL, EVAL1)
    }
  }
  
  colnames(DATA) <- gsub(paste("RUN", NbRunEval + 1, sep = ""), 
                         "Full", colnames(DATA))
  colnames(DATA)[3:ncol(DATA)] <- paste(colnames(DATA)[3:ncol(DATA)], 
                                        "ESM", sep = "_")
  DATA[, -c(1, 2)] <- DATA[, -c(1, 2)] * 1000
  EVAL[, 2:14] <- round(EVAL[, 2:14], 3)
  EVAL[, 2] <- EVAL[, 2] * 1000
  if (length(models) > 1) {
    output <- list(species = data@sp.name, ESM.fit = round(DATA[, 
                                                                -1]), ESM.evaluations = EVAL, weights = weights, 
                   weights.EF = weights.double, failed = failed.mod)
  }
  else {
    output <- list(species = data@sp.name, ESM.fit = round(DATA[, 
                                                                -1]), ESM.evaluations = EVAL, weights = weights, 
                   failed = failed.mod)
  }
  save(output, file = paste("ESM_EnsembleModeling..", weighting.score, 
                            threshold, ESM.modeling.output$modeling.id, "out", sep = "."))
  setwd(iniwd)
  return(output)
}

############################
# d)  ecospat.ESM.EnsembleProjection
############################


#Description
## Projecting calibrated ESMs into new space or time

## FUNCTION'S ARGUMENTS
## ESM.prediction.output,:          Object returned by ecospat.ESM.Projection
## ESM.EnsembleModeling.output:     Object returned by ecospat.ESM.EnsembleModeling

## Values:
# ESM.projections:    Returns the projections of ESMs for the selected single models and their ensemble (data frame or SpatRaster).


## Authors:
# Frank Breiner
# With the updates of Flavien Collart

## References
#Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
# Breiner F.T., Guisan A., Bergamini A., Nobis M.P. (2015):
#  Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218. doi: 10.1111/2041-210X.12403

## See also
# ecospat.ESM.Modeling; ecospat.ESM.MergeModels; ecospat.ESM.EnsembleModeling

ecospat.ESM.EnsembleProjection <- function(ESM.prediction.output, ESM.EnsembleModeling.output, chosen.models = 'all') {
  models <- ESM.prediction.output$models
  weights <- ESM.EnsembleModeling.output$weights
  weights.EF <- ESM.EnsembleModeling.output$weights.EF
  if (chosen.models[1] != "all") {
    if (any(!chosen.models %in% models)) {
      stop("chosen.models must be a subset of the models selected in ecospat.ESM.Modeling()")
    }
    models <- chosen.models
    weights.EF <- weights.EF[weights.EF$Group.1 %in% chosen.models, 
    ]
  }
  NbRunEval <- ESM.prediction.output$NbRunEval
  pred.biva <- ESM.prediction.output$pred.biva
  new.env.raster <- ESM.prediction.output$new.env.raster
  failed.mod <- grep("allRun", unlist(ESM.EnsembleModeling.output$failed),
                     value = TRUE)
  if (length(models) == 1) {
    weigths.rm <- NULL
    for (i in 1:length(weights)) {
      weigths.rm <- c(weigths.rm, which(grepl(paste(names(weights), 
                                                    ".", sep = "")[i], pred.biva, fixed = TRUE)))
    }
    pred.biva <- pred.biva[weigths.rm]
    pred.biva <- sort(pred.biva)
  }
  if (new.env.raster) {
    pred.biva <- grep("\\.tif\\b", pred.biva, value = TRUE)
    pred.biva <- pred.biva[gtools::mixedorder(gsub(".", "_", pred.biva, 
                                           fixed = TRUE))]
  }
  if (!new.env.raster) {
    string_pred <- grep("RData", pred.biva, value = TRUE)
    pred.biva <- string_pred[gtools::mixedorder(gsub(".", "_", string_pred, 
                                             fixed = TRUE))]
  }
  biva.proj <- list()
  for (i in 1:length(pred.biva)) {
    if (!new.env.raster) {
      
      aa <- as.data.frame(get(load(pred.biva[i])))
      nameMap <- gsub(".*proj_","",pred.biva[i])
      nameMap <- sub(paste0(ESM.prediction.output$name.env,"."),"",nameMap,fixed = T)
      nameMap <- sub(paste0(".",ESM.prediction.output$modeling.id,".RData"),"",nameMap,fixed = T)
      colnames(aa) = nameMap
      biva.proj[[i]] <- aa
      
    }
    if (new.env.raster) {
      biva.proj[[i]] <- terra::rast(pred.biva[i])
      nameMap <- gsub(".*proj_","",pred.biva[i])
      nameMap <- sub(paste0(ESM.prediction.output$name.env,"."),"",nameMap,fixed = T)
      nameMap <- sub(paste0(".",ESM.prediction.output$modeling.id,".tif"),"",nameMap,fixed = T)
      names(biva.proj[[i]]) = nameMap
    }
  }
  if (!new.env.raster) {
    pred.ESM <- list()
    biva.st <- do.call(cbind, biva.proj)
    for (i in 1:length(models)) {
      Model.biva.prediction <- biva.st[, grep(models[i], 
                                              colnames(biva.st))]
      
      weights2 <- weights[colnames(Model.biva.prediction)]
      pred.ESM[[i]] <- apply(Model.biva.prediction, 1, 
                             function(x) stats::weighted.mean(x, weights2, 
                                                              na.rm = TRUE))
    }
    rm(biva.st)
    pred.ESM <- round(as.data.frame(do.call(cbind, pred.ESM)))
    colnames(pred.ESM) <- models
  }
  # if (length(models) == 1) {
  #   names(weights) <- paste0(models, names(weights))
  # }
  if (new.env.raster) {
    pred.ESM <- terra::rast(biva.proj)
    for (i in 1:length(models)) {
      model.maps<- pred.ESM[[grep(models[i], names(pred.ESM))]]
      weights.mod <- weights[grep(models[i], names(weights))]
      if (any(grepl(models[i], failed.mod))) {
        weights.mod <- weights.mod[!names(weights.mod) %in% 
                                     paste(models[i], ".", sub("_.*", "", failed.mod[grepl(models[i], 
                                                                                           failed.mod)]), sep = "")]
      }
      model.maps <- model.maps[[weights.mod!=0]] 
      weights.mod <- weights.mod[weights.mod>0]
      assign(models[i], round(terra::weighted.mean(model.maps,
                                                   weights.mod, na.rm = TRUE)))
    }
    pred.ESM <- terra::rast(mget(models))[[order(models)]]
    do.call("rm", as.list(models))
  }
  if (length(models) > 1) {
    if (new.env.raster) {
      ESM.EF <- round(terra::weighted.mean(pred.ESM[[names(pred.ESM)]], 
                                            weights.EF[order(weights.EF[, 1]), 2], na.rm = TRUE))
      pred.ESM <- c(pred.ESM, ESM.EF)
      names(pred.ESM) <- c(names(pred.ESM)[1:(terra::nlyr(pred.ESM) - 
                                                1)], "EF")
      rm(ESM.EF)
    }
    if (!new.env.raster) {
      pred.ESM$EF <- apply(pred.ESM, 1, function(x) stats::weighted.mean(x, 
                                                                         weights.EF[order(weights.EF[, 1]), 2], na.rm = TRUE))
    }
    if (!new.env.raster) {
      pred.ESM <- round(pred.ESM)
    }
  }
  return(ESM.projections = pred.ESM)
}
############################
# e)  ecospat.mpa
############################

## This function calculates the Minimal Predicted Area.

## FUNCTION'S ARGUMENTS
## Pred:      numeric or, SpatRaster .predicted suitabilities from a SDM prediction
## Sp.occ.xy: xy-coordinates of the species (if Pred is a  SpatRaster)
## perc:      Percentage of Sp.occ.xy that should be encompassed by the binary map.

## Details:

## Value
# Returns the Minimal Predicted Area

## Author(s)
# Frank Breiner with the contribution of Flavien Collart

## References
# Engler, R., A. Guisan, and L. Rechsteiner. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. Journal of Applied Ecology 41:263-274.


ecospat.mpa <- function(Pred, Sp.occ.xy = NULL, perc = 0.9) {
  perc <- 1 - perc
  if (!is.null(Sp.occ.xy)) {
    Pred <- terra::extract(Pred, as.data.frame(Sp.occ.xy),ID=FALSE)
  }
  round(quantile(Pred, probs = perc,na.rm = TRUE), 3)
}

### EXAMPLE

# obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
# [which(ecospat.testData$Saxifraga_oppositifolia==1)]) ecospat.mpa(obs) ecospat.mpa(obs,perc=1) ##
# 100% of the presences encompassed


## Example script for using Ensemble Small Models ESMs according to Lomba et al. 2010 Written by
## Frank Breiner Swiss Federal Research Institute WSL, July 2014.

## References: Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J. &
## Guisan, A. (2010). Overcoming the rare species modeling paradox: A novel hierarchical framework
## applied to an Iberian endemic plant. Biological Conservation, 143, 2647-2657.



### f) ecospat.ESM.MergeModels merge different ecospat.ESM.Modeling outputs


# Description Enables to merge different ecospat.ESM.Modeling outputs which were produced e.g. on
# different computers or on a cluster


## FUNCTION'S ARGUMENTS ESM.modeling.output: list object of different ESM.modeling.output objects.

# Details The Biomod output folders must be copied into one working directory if the output is at
# different locations!


## Values: modeling.id: 'character', id of modeling procedure models.: location of models output
## models: models argument ESM_Projection: logical. TRUE if models were projected for new.env
## pred.biva: location of the bivariate models predictions calib.lines: 'BIOMOD.stored.array',
## calibrations and evaluation lines selected for each run NbRunEval: Number of Evaluation run
## new.env: set of explanatory variables onto which models were projected data:
## 'BIOMOD.formated.data' object wd: working directory which.biva: Which bivariate combinations were
## used for modeling

## Note: You should use this function only if you want to split your bivariate model runs, e.g. on a
## clusters.

## See also: ecospat.ESM.Modeling, ecospat.ESM.EnsembleModeling

## Authors: Frank Breiner


ecospat.ESM.MergeModels <- function(ESM.modeling.output) {
  number.outputs <- length(ESM.modeling.output)
  
  ## tests
  if (length(unique(sapply(ESM.modeling.output, function(x) x[3]))) != 1) {
    stop("models of ESM.modeling.output differ")
  }
  if (length(unique(sapply(ESM.modeling.output, function(x) x[7]))) != 1) {
    stop("NbRunEval of ESM.modeling.output differ")
  }
  if (length(unique(sapply(ESM.modeling.output, function(x) x[8]))) != 1) {
    stop("new.env of ESM.modeling.output differ")
  }
  if (length(unique(sapply(ESM.modeling.output, function(x) x[9]))) != 1) {
    stop("data of ESM.modeling.output differ")
  }
  # if(length(unique(sapply( ESM.modeling.output,function(x) x[2])))==1){stop('The ESM models you want
  # to merge are all identical')}
  combinations <- combn(colnames(ESM.modeling.output[[1]]$data@data.env.var), 2)
  test <- unique(do.call(c, sapply(ESM.modeling.output, function(x) x[11])), MARGIN = 2)
  # miss <-
  # !duplicated(cbind(test,combinations),MARGIN=2)[(ncol(test)+1):(ncol(test)+ncol(combinations))]
  if (ncol(combinations) != length(test)) {
    warning(cat(paste("\n\n\n################################\n", "There are", ncol(combinations) -
                        length(test), "Biva models missing for the ESM. You could run ecospat_ESM_Modeling() with which.biva = ",
                      list(which(!1:ncol(combinations) %in% test)), " to include them.", "\n################################\n",
                      sep = " ")))
    warning("missing: ", print(combinations[, which(!1:ncol(combinations) %in% test)]))
  }
  ESM.modeling.output[[1]]$modeling.id <- do.call(c, unique(sapply(ESM.modeling.output, function(x) x[1])))
  ESM.modeling.output[[1]]$models. <- unique(do.call(c, unique(sapply(ESM.modeling.output, function(x) x[2]))))
  ESM.modeling.output[[1]]$pred.biva <- unique(do.call(c, unique(sapply(ESM.modeling.output, function(x) x[5]))))
  ESM.modeling.output[[1]]$which.biva <- unique(do.call(c, unique(sapply(ESM.modeling.output, function(x) x[11]))))
  
  del.pred <- select.pred <- del.models <- select.models <- NULL
  for (i in 1:length(ESM.modeling.output[[1]]$which.biva)) {
    select.models <- c(select.models, grep(paste("ESM.BIOMOD.", ESM.modeling.output[[1]]$which.biva[i],
                                                 ".", sep = ""), ESM.modeling.output[[1]]$models., fixed = TRUE)[1])
    del.models <- c(del.models, grep(paste("ESM.BIOMOD.", ESM.modeling.output[[1]]$which.biva[i],
                                           ".", sep = ""), ESM.modeling.output[[1]]$models., fixed = TRUE)[-1])
    
    if (length(grep(".tif", ESM.modeling.output[[1]]$pred.biva)) > 0) {
      select.pred <- c(select.pred, grep(paste("ESM.BIOMOD.", ESM.modeling.output[[1]]$which.biva[i],
                                               ".tif", sep = ""), ESM.modeling.output[[1]]$pred.biva)[1])
      del.pred <- c(del.pred, grep(paste("ESM.BIOMOD.", ESM.modeling.output[[1]]$which.biva[i],
                                         ".tif", sep = ""), ESM.modeling.output[[1]]$pred.biva)[-1])
    }
    if (length(grep(".RData", ESM.modeling.output[[1]]$pred.biva)) > 0) {
      select.pred <- c(select.pred, grep(paste("ESM.BIOMOD.", ESM.modeling.output[[1]]$which.biva[i],
                                               ".RData", sep = ""), ESM.modeling.output[[1]]$pred.biva)[1])
      del.pred <- c(del.pred, grep(paste("ESM.BIOMOD.", ESM.modeling.output[[1]]$which.biva[i],
                                         ".RData", sep = ""), ESM.modeling.output[[1]]$pred.biva)[-1])
    }
  }
  if (length(del.models) > 0) {
    warning(cat(paste("\n\n\n################################\n", length(del.models), "models are not considered for ESMs because they seem to be either replicates or artefacts (e.g. BIOMOD.models.out-objects which already existed before running ecospat_ESM_Modeling function but do not belong to the ESM).\n
                      The following models were deleted:\n",
                      sep = " ")), print(ESM.modeling.output[[1]]$models.[del.models]))
  }
  if (length(del.pred) > 0) {
    warning(cat(paste("\n\n\n################################\n", length(del.pred), "predictions were deleted because they seem to be replicates:\n",
                      sep = " ")), print(ESM.modeling.output[[1]]$pred.biva[del.pred]))
  }
  ESM.modeling.output[[1]]$models. <- ESM.modeling.output[[1]]$models.[select.models]
  ESM.modeling.output[[1]]$pred.biva <- ESM.modeling.output[[1]]$pred.biva[select.pred]
  
  return(ESM.modeling.output[[1]])
}

############################
# b)  ecospat.ESM.VarContrib          get Variable contribution for each variable and method (mean across model runs)
############################

## FUNCTION'S ARGUMENTS
## ESM.modeling.output:   BIOMOD.formated.data object returned by BIOMOD_modeling
## ESM_EF.output:   BIOMOD.formated.data object returned by ecospat.ESM.EnsembleModeling

## Details:
# Calculates the ratio between sum of weights of bivariate models where a focal variable was used and sum of weights of all bivariate models. This gives an indication on the proportional contribution of the variable in the final ensemble model. 
# In the case of multiple methods (e.g., GLM, GAM...), the contributions are counted per method, and for the final ensemble model, the contributions are then weighted means of single methods (based on the weighting score as chosen in ecospat.ESM.EnsembleModeling()).

## Values:
# Returns a dataframe with contribution values (i.e., proportional contribution) by variable and model

## Authors:
# Olivier Broennimann <Olivier.Broennimann@unil.ch>

##See Also
# ecospat.ESM.Modeling; ecospat.ESM.EnsembleModeling; ecospat.ESM.EnsembleProjection

ecospat.ESM.VarContrib <- function(ESM.modeling.output,ESM_EF.output) {
  
  var<-colnames(ESM.modeling.output$data@data.env.var)
  models<-ESM.modeling.output$models
  contrib<-data.frame(matrix(nrow=length(var),ncol=length(models),dimnames=list(var,c(models))))
  weights<-ESM_EF.output$weights
  if(length(models)==1){
    names(weights) <- paste0(models,names(weights))
    order_weights <- order(as.numeric(sub("GLM.ESM.BIOMOD.","",names(weights))))
  }
  
  if(length(models)>1){
    order_weights <- paste0(rep(models, ncol(combn(var,2))), ".ESM.BIOMOD.", 
                            rep(1:ncol(combn(var,2)), each=length(models)))
  }
  weights.reordered<-weights[order_weights]
  
  #contribution by technique
  
  cb1<-rep(combn(var,2)[1,],each=length(models))
  cb2<-rep(combn(var,2)[2,],each=length(models))
  
  for (m in models){
    for(v in var){
      pos_models <- grep(m,names(weights.reordered))
      pos_cb<-c(which(cb1==v),which(cb2==v))
      pos<-intersect(pos_models,pos_cb)
      contrib[which(var==v),which(names(contrib)==m)] <- 
        sum(weights.reordered[pos])/(sum(weights.reordered[setdiff(pos_models,pos)]))*length(setdiff(pos_models,pos))/length(pos)
      #<-sum(weights.reordered[pos])/(2*sum(weights.reordered[pos_models]))
    }
  }
  
  #contributions of final ensemble model
  if(length(models) > 1) {
    EF <- matrixStats::rowWeightedMeans(x=data.matrix(contrib[, models]), w=ESM_EF.output$weights.EF$x, na.rm=TRUE)
    contrib<-cbind(contrib,EF) }
  
  return(contrib)
}

 #Function to generate the argument DataSplitTable in the function ecospat.ESM.Modeling                                                                        
.CreatingDataSplitTable <- function(bm.format, 
                                    NbRunEval, 
                                    DataSplit){
  xy <- bm.format@coord
  resp <- bm.format@data.species
  pres <- which(resp==1 & !is.na(resp))
  abs <- setdiff(1:length(resp),pres)
  calib.Lines <- matrix(FALSE, nrow = length(resp), ncol = NbRunEval)
  if ("PA.table" %in% slotNames(bm.format)) {
    abs <- intersect(abs,which(bm.format@PA.table[,1]))
  }
  
  for(i in 1:NbRunEval){
    calib.Lines[sample(pres,size = round(length(pres)*DataSplit/100)),i] = TRUE
    calib.Lines[sample(abs,size = round(length(abs)*DataSplit/100)),i] = TRUE
  }
  calib.Lines <- cbind(calib.Lines,TRUE)
  colnames(calib.Lines) = c(paste0("_RUN",1:NbRunEval),"_Full")
  if ("PA.table" %in% slotNames(bm.format)) {
    calib.Lines[!(bm.format@PA.table[,1]),] = NA
  }
  return(calib.Lines)
}
