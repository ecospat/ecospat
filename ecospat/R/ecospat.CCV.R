#####################################################################
# Function for community cross-validation modelling and evaluation
#####################################################################
# Commented version of the ecospat.CCV functions to model and evaluate 
# species assemblages (Scherrer et al. 2018) using normals species distribution models (SDMs; biomod2)
# or an ensemble of small models (ESMs; ecospat.ESM).

# The ecospat.CCV group contains four main function to

# a) create a cross-validation dataset to be used (ecospat.CCV.createDataSplitTable):
#    This function creates a DataSplitTable with calibration and evaluation data
#    either for cross-validation or repeated split sampling at the community level
#    (i.e., across all species)
#
# b) run the models for all individual species of the community (ecospat.CCV.modeling)
#    This function runs individual SDMs () and/or ESMs (Breiner et al. 2015) for all the species in a community
#    and returns their predictions (i.e., probabilities/habitat suitability) along with
#    a range of evaluation metrics for the models.
#
# c) evaluate the community predictions based on different thresholding techniques (ecospat.CCV.communityEvaluation.bin)
#    This function binarised the individual predictions and then stacks the binary maps. There are a range
#    of different binarisation methods available and the resulting species assemblages get evaluated 
#    by different community evaluation metrics
#
# d) evaluate the community predictions based directly on the probabilites (i.e., threshold indpendent) (ecospat.CCV.communityEvaluation.prob)
#    This function generates a number of community evaluation metrics directly based on the probability/habitat suitability
#    returned by the individual models.
#
#
#
# REFERENCES:
# Scherrer, D., D'Amen, M., Mateo, M.R.G., Fernandes, R.F. & Guisan , A. (2018) 
#   How to best threshold and validate stacked species assemblages? Community optimisation might hold the answer. 
#   Methods in Ecology and Evolution, 9(10): 2155-2166.
# Thuiller, W., Georges, D., Engler, R. & Breiner, F. (2016) 
#   biomod2: Ensemble Platform for Species Distribution Modeling.
# Breiner, F.T., Guisan, A., Bergamini, A. & Nobis, M.P. (2015) 
#   Overcoming limitations of modelling rare species by using ensembles of small models. 
#   Methods in Ecology and Evolution, 6, 1210-1218.
#
#
# Daniel Scherrer, University of Lausanne, daniel.j.a.scherrer@gmail.com, July 2018
##################################################################################################################




#####################################################
# a) ecospat.CCV.createDataSplitTable
#####################################################

#DESCRIPTION
# Creates a DataSplitTable with calibration and evaluation data
# either for cross-validation or repeated split sampling at the community level
# (i.e., across all species)


#FUNCTION'S ARGUMENTS
#NbSites              number of total sites available
#NbRunEval            number of cross-validation or split sample runs 
#DataSplit            proporation of sites used for model calibration 
#validation.method    the type of DataSplitTable that should be created. Must be either "cross-validation" or "split-sample".


#DETAILS


#VALUES
#DataSplitTable       a matrix with TRUE/FALSE for each model run (TRUE=Calibration point, FALSE=Evaluation point)


#AUTHORS
# Daniel Scherrer <daniel.j.a.scherrer@gmail.com>


ecospat.CCV.createDataSplitTable <- function(NbRunEval, 
                                             DataSplit,
                                             validation.method,
                                             NbSites,
                                             sp.data=NULL,
                                             minNbPresences=NULL,
                                             minNbAbsences=NULL,
                                             maxNbTry=1000){
  
  #Check the input data
  stopifnot(DataSplit >= 50 & DataSplit <=100)
  stopifnot(NbRunEval>=0)
  stopifnot(validation.method %in% c("cross-validation", "split-sample"))
  if(is.null(sp.data)){
    stopifnot(NbSites>0)
  }else{
    stopifnot(!is.null(minNbPresences) & !is.null(minNbAbsences) & minNbPresences >= 0 & minNbAbsences >= 0 & maxNbTry>0 & maxNbTry < 1000000000)
    stopifnot(is.data.frame(sp.data))
  }
  
  
  
  #################################################
  ## Simple version without species data check ####
  #################################################
  
  if(is.null(sp.data)){
    #Create the empty table
    DataSplitTable <- matrix(data=FALSE, nrow=NbSites, ncol=NbRunEval)
    
    #Split-sample appraoch
    if(validation.method=="split-sample"){
      for(i in 1:NbRunEval){
        DataSplitTable[sample(1:NbSites,round(DataSplit/100*NbSites), replace=FALSE),i] <- TRUE
      }
    }
    #Cross-validation approach
    if(validation.method=="cross-validation"){
      grouper <- sample(rep(1:NbRunEval,each=ceiling(NbSites/NbRunEval)),NbSites, replace = FALSE)
      iner <- round(DataSplit*NbRunEval/100)
      for(i in 1:NbRunEval){
        DataSplitTable[which(grouper %in% ((i:(i+iner-1)%%NbRunEval)+1)),i] <- TRUE
      }
    }
    return(DataSplitTable)
  }
  
  #################################################
  ## Complex version with species data check     ##
  #################################################
  
  if(!is.null(sp.data)){
    
    #Create Species x Run matrix ################################################################################
    create.SpRunMatrix <- function(sp.data, 
                                   DataSplitTable){
      SpRunMatrix <- apply(DataSplitTable,2, function(x){colSums(x*sp.data)})
    }
    
    Nb.sp.dropped <- dim(sp.data)[2]
    trys <- 1
    
    while(trys <= maxNbTry & Nb.sp.dropped > 0){
      #Create the empty table
      DataSplitTable <- matrix(data=FALSE, nrow=dim(sp.data)[1], ncol=NbRunEval)
      
      #Split-sample appraoch
      if(validation.method=="split-sample"){
        for(i in 1:NbRunEval){
          DataSplitTable[sample(1:dim(sp.data)[1],round(DataSplit/100*dim(sp.data)[1]), replace=FALSE),i] <- TRUE
        }
      }
      
      #Cross-validation approach
      if(validation.method=="cross-validation"){
        grouper <- sample(rep(1:NbRunEval,each=ceiling(dim(sp.data)[1]/NbRunEval)),dim(sp.data)[1], replace = FALSE)
        iner <- round(DataSplit*NbRunEval/100)
        for(i in 1:NbRunEval){
          DataSplitTable[which(grouper %in% ((i:(i+iner-1)%%NbRunEval)+1)),i] <- TRUE
        }
      }
      
      #SpRunMatrix
      SpRunMatrix <- create.SpRunMatrix(sp.data = sp.data, DataSplitTable = DataSplitTable)
      
      #Check for the species
      sp.names.all <- rownames(SpRunMatrix)
      sp.names.ok <- intersect(names(which(apply(SpRunMatrix,1,min) >= minNbPresences)), names(which(apply(min(colSums(DataSplitTable))-SpRunMatrix,1,min) >= minNbAbsences)))
      sp.names.droped <- setdiff(sp.names.all, sp.names.ok)
      
      #Save the best DataSplitTable
      if(length(sp.names.droped) < Nb.sp.dropped){
        DataSplitTable.final <- DataSplitTable
        Nb.sp.dropped <- length(sp.names.droped)
        sp.names.droped.best <- sp.names.droped
      }
      
      trys <- trys+1
    }
    message(paste("The following species will not have the desired minimum number of presence/absence data in each run: ", paste(sp.names.droped.best, sep="", collapse=", "),"\n\n",sep=""))
    return(DataSplitTable.final)
  }
}

#####################################################
# b) ecospat.CCV.modeling
#####################################################

#DESCRIPTION
# Creates probabilistic prediction for all species based on SDMs or ESMs 
# and returns their evaluation metrics and variable importances.


#FUNCTION'S ARGUMENTS
#sp.data              a binary matrix where the rows are sites and the columns are species (values 1,0)
#env.data             either a data.frame where rows are sites and colums are environmental variables or a raster stack of the envrionmental variables
#xy                   2 comumn data.frame with X and Y coordinates of the sites (most be same coordinate system as env.data)
#DataSplitTable       a table providing TRUE/FALSE to indicate what points are used for calibration and evaluation. As returned by ecospat.CCV.createDataSplitTable().
#DataSplit            proporation of the data used for model calibration (only needed if no DataSplitTable provided)
#NbRunEval            number of cross-validatio/split sample runs (only needed if no DataSplitTable provided)
#minNbPredictors      minimum number of occurences [min(presences/Absences] per predicotors needed to calibrate the models
#validation.method    either "cross-validation" or "split-sample" used to validate the communtiy predictions (only needed if no DataSplitTable provided)
#models.sdm           modeling techniques used for the normal SDMs. Vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips' and 'MAXENT.Tsuruoka'
#models.esm           modeling techniques used for the ESMs
#modeling.options.sdm BIOMOD.models.options object returned by BIOMOD_ModelingOptions (same as in biomod2) for normal SDMs
#modeling.options.esm BIOMOD.models.options object returned by BIOMOD_ModelingOptions (same as in biomod2) for esm SDMs
#ensemble.metric      metrics used to create the ensemble models
#ESM                  either YES (ESMs allowed), NO (ESMs not allowed) or "ALL" (ESMs used in any case)
#parallel             should parallel computing be allowed
#cpus                 number of cpus to use in parallel computing
#VarImport            number of permutation runs to evaluate variable importance
#modeling.id          name of the model


#DETAILS


#VALUES
#modelling.id                                             "character string" of the model name
#output.files                                             vector with the file names written to the hard drive
#speciesData.calibration                                  a 3-dimensional array of presence/absence data of all species for the calibration plots used for each run
#speciesData.evaluation                                   a 3-dimensional array of presence/absence data of all species for the evaluation plots used for each run
#speciesData.full                                         a matrix of presence/absence data of all species
#DataSplitTable                                           a matrix with TRUE/FALSE for each model run (TRUE=Calibration point, FALSE=Evaluation point)
#singleSpecies.ensembleEvaluationScore                    a 3-dimensional array of single species evaluation metrics (Max.KAPPA, Max.TSS, AUC of ROC)
#singleSpecies.ensembleVariableImportance                 a 3-dimensional array of single species variable importance for all predictors
#singleSpecies.calibrationSites.ensemblePredictions       a 3-dimensional array of the predictions for each species and run at the calibration sites
#singleSpecies.evaluationSites.ensemblePredictions        a 3-dimensional array of the predictions for each species and run at the evaluation sites


#AUTHORS
# Daniel Scherrer <daniel.j.a.scherrer@gmail.com>


#REFERENCES
# Scherrer, D., D'Amen, M., Mateo, M.R.G., Fernandes, R.F. & Guisan , A. (2018) 
#   How to best threshold and validate stacked species assemblages? Community optimisation might hold the answer. 
#   Methods in Ecology and Evolution, 9(10): 2155-2166.
# Thuiller, W., Georges, D., Engler, R. & Breiner, F. (2016) 
#   biomod2: Ensemble Platform for Species Distribution Modeling.
# Breiner, F.T., Guisan, A., Bergamini, A. & Nobis, M.P. (2015) 
#   Overcoming limitations of modelling rare species by using ensembles of small models. 
#   Methods in Ecology and Evolution, 6, 1210-1218.


#SEE ALSO
# ecospat.CCV.createDataSplitTable, ecospat.CCV.communityEvaluation.bin, ecospat.CCV.communityEvaluation.prob


ecospat.CCV.modeling <- function(sp.data, 
                                 env.data, 
                                 xy,
                                 DataSplitTable=NULL,
                                 DataSplit = 70, 
                                 NbRunEval = 25,
                                 minNbPredictors =5,
                                 validation.method = "cross-validation",
                                 models.sdm = c("GLM","RF"), 
                                 models.esm = "CTA", 
                                 modeling.options.sdm = NULL, 
                                 modeling.options.esm = NULL, 
                                 ensemble.metric = "AUC", 
                                 ESM = "YES",
                                 parallel = FALSE, 
                                 cpus = 4,
                                 VarImport = 10,
                                 modeling.id = as.character(format(Sys.time(), '%s'))){
  
  #Loading the packages needed (they should all be installed by ecospat library)
  #require(raster)
  #require(rgdal)
  #require(biomod2)
  #require(snowfall)
  #require(gtools)
  #require(PresenceAbsence)
  
  #Checking all the input data
  stopifnot(dim(sp.data)[1]==dim(xy)[1])
  stopifnot(dim(env.data)[1]==dim(xy)[1] | data.class(env.data)=="RasterStack")
  stopifnot(dim(DataSplitTable)[1]==dim(xy)[1] | is.null(DataSplitTable))
  stopifnot(DataSplit >= 50 & DataSplit <=100)
  stopifnot(NbRunEval>=0)
  stopifnot(minNbPredictors>1)
  stopifnot(validation.method %in% c("cross-validation", "split-sample"))
  stopifnot(models.sdm %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips','MAXENT.Tsuruoka'))
  stopifnot(models.esm %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips','MAXENT.Tsuruoka'))
  stopifnot(ensemble.metric %in% c("AUC","TSS","KAPPA") & length(ensemble.metric)==1)
  stopifnot(ESM %in% c("YES","NO","ALL"))
  stopifnot(is.logical(parallel))
  stopifnot(cpus>=1)
  stopifnot(is.numeric(VarImport))
  
  #Create the variables for individual model evaluation and assembly
  eval.metrics.sdm=c('KAPPA', 'TSS', 'ROC')
  eval.metrics.esm=c('KAPPA', 'TSS', 'AUC')
  eval.metrics.names= c('KAPPA', 'TSS', 'AUC')
  ensemble.metric.esm <- ensemble.metric
  if(ensemble.metric=="AUC"){
    ensemble.metric.sdm <- "ROC"
  }else{
    ensemble.metric.sdm <- ensemble.metric
  }
  if(data.class(env.data)=="RasterStack"){
    NbPredictors <- dim(env.data)[3]
    NamesPredictors <- names(env.data)
  }else{
    NbPredictors <- dim(env.data)[2]
    NamesPredictors <- colnames(env.data)
  }
  if(length(models.esm)==1){
    ef.counter <- 1
  }else{
    ef.counter <- length(models.esm)+1
  }
  
  #Adjust species names to have dots
  colnames(sp.data) <- gsub("_",".", colnames(sp.data))
  
  #Create the folder to write the data output
  dir.create(modeling.id)
  oldwd <- getwd() 
  on.exit(setwd(oldwd)) 
  setwd(modeling.id)
  
  #############################################################################################################
  #SUBFUNCTIONS################################################################################################
  #############################################################################################################
  
  #Create Species x Run matrix ################################################################################
  create.SpRunMatrix <- function(sp.data, 
                                 DataSplitTable){
    SpRunMatrix <- apply(DataSplitTable,2, function(x){colSums(x*sp.data)})
  }
  
  #Function to run biomod2 in parallel ########################################################################
  BiomodSF <- function(sp.name, 
                       DataSplitTable, 
                       sp.data, 
                       env.data, 
                       xy, 
                       models, 
                       models.options, 
                       eval.metrics, 
                       ensemble.metric,
                       VarImport){
    
    #Preparing the data
    MyBiomodData <- BIOMOD_FormatingData(resp.var = as.numeric(sp.data[,sp.name]),
                                         expl.var = env.data,
                                         resp.xy = xy,
                                         resp.name = sp.name,
                                         na.rm=FALSE)
    
    #Setting model parameters
    if(is.null(models.options)){
      MyBiomodOptions <- BIOMOD_ModelingOptions()
    }else{
      MyBiomodOptions <- BIOMOD_ModelingOptions(models.options) #NOT WORKING YET!!!!
    }
    
    #Running the models
    MyBiomodModelOut <- BIOMOD_Modeling(data = MyBiomodData,
                                        models = models,
                                        models.options = MyBiomodOptions,
                                        models.eval.meth = eval.metrics,
                                        DataSplitTable = DataSplitTable,
                                        Prevalence=NULL,
                                        modeling.id = "ccv")
    
    #Creating the ensemble Model
    MyBiomodEnsemble <- BIOMOD_EnsembleModeling(modeling.output = MyBiomodModelOut,
                                                chosen.models = "all",
                                                em.by = "PA_dataset+repet",
                                                eval.metric = ensemble.metric,
                                                eval.metric.quality.threshold = NULL,
                                                models.eval.meth =eval.metrics,
                                                prob.mean = FALSE,
                                                prob.cv = FALSE,
                                                prob.ci = FALSE,
                                                prob.ci.alpha = 0.05,
                                                prob.median = FALSE,
                                                committee.averaging = FALSE,
                                                prob.mean.weight = TRUE,
                                                prob.mean.weight.decay = 'proportional',
                                                VarImport = VarImport)
  }
  
  #Function to run ESM in parallel #######################################################################
  ESMSF <- function(sp.name, 
                    DataSplitTable, 
                    sp.data, 
                    env.data, 
                    xy, 
                    models, 
                    models.options, 
                    ensemble.metric){
    
    #Preparing the data
    MyESMData <- BIOMOD_FormatingData(resp.var = as.numeric(sp.data[,sp.name]),
                                      expl.var = env.data,
                                      resp.xy = xy,
                                      resp.name = sp.name,
                                      na.rm = FALSE)
    
    #Setting model parameters
    if(is.null(models.options)){
      MyBiomodOptions <- BIOMOD_ModelingOptions()
    }else{
      MyBiomodOptions <- BIOMOD_ModelingOptions(models.options) #NOT WORKING YET!!!!
    }
    
    #Running the ESM
    MyESMModelOut <- ecospat.ESM.Modeling(data=MyESMData, 
                                          DataSplitTable = DataSplitTable, 
                                          weighting.score = ensemble.metric,
                                          models=models,
                                          Prevalence=NULL,
                                          modeling.id="ccv", 
                                          models.options=MyBiomodOptions, 
                                          parallel=FALSE)
    
    #Ensemble the ESMs
    MyESMEnsemble <- ecospat.ESM.EnsembleModeling(ESM.modeling.output = MyESMModelOut,
                                                  weighting.score = ensemble.metric,
                                                  models=models)
  }
  
  #Get the variable contribution form the ESM models
  get.ESMvariableContribution <- function(output_EF, output, NamesPredictors){
    
    #Dataframe to save the variable contribution
    Variable.Contribution <- rep(NA, times=length(NamesPredictors))
    names(Variable.Contribution)<- NamesPredictors
    
    #Calculating the contribution for each predictor
    for(v in NamesPredictors){
      cb1<-rep(combn(NamesPredictors,2)[1,],each=length(output$models))
      cb2<-rep(combn(NamesPredictors,2)[2,],each=length(output$models))
      pos<-c(which(cb1==v),which(cb2==v))
      Variable.Contribution[which(NamesPredictors==v)]<-mean(output_EF$weights[pos])-mean(output_EF$weights)
    }
    
    #Returning the output
    Variable.Contribution[which(is.na(Variable.Contribution))]<-0 #to remove variables with no contribution NA
    return((Variable.Contribution-min(Variable.Contribution))/(max(Variable.Contribution)-min(Variable.Contribution)))
  }
  
  ##############################################################################################################
  ##############################################################################################################
  ##############################################################################################################
  
  #Creating the calibration and evaluation data (if not provided)
  if(is.null(DataSplitTable)){
    DataSplitTable <- ecospat.CCV.createDataSplitTable(NbSites = dim(xy)[1], NbRunEval = NbRunEval, DataSplit = DataSplit, validation.method = validation.method)
  }else{
    NbRunEval <- dim(DataSplitTable)[2]
  }
  
  #Creating a speciesxrun matrix
  SpRunMatrix <- create.SpRunMatrix(sp.data = sp.data, DataSplitTable = DataSplitTable)
  
  #Running the SDMs for the species
  
  #Running the stuff with NO ESMs###############################################################################
  if(ESM=="NO"){
    
    #Selecting the species than can be run (enough data)
    sp.names.all <- rownames(SpRunMatrix)
    sp.names.ok <- intersect(names(which(apply(SpRunMatrix,1,min) >= minNbPredictors*NbPredictors)), names(which(apply(min(colSums(DataSplitTable))-SpRunMatrix,1,min) >= minNbPredictors*NbPredictors)))
    sp.names.droped <- setdiff(sp.names.all, sp.names.ok) 
    message(paste("The following species will not be modelled due to limited presence data: ", paste(sp.names.droped, sep="", collapse=", "),"\n\n",sep=""))
    
    #Creating the speciesData.calibration and speciesData.evaluation
    speciesData.calibration <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(DataSplitTable)),  NbRunEval), dimnames=list(sp.names.ok, paste("c",1:max(colSums(DataSplitTable)),sep="_"), 1:NbRunEval))
    speciesData.evaluation <- singleSpecies.evaluationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(!DataSplitTable)), NbRunEval), dimnames=list(sp.names.ok, paste("e",1:max(colSums(!DataSplitTable)),sep="_"), 1:NbRunEval))
    for(i in 1:NbRunEval){
      speciesData.calibration[,1:sum(DataSplitTable[,i]),i] <- t(sp.data[which(DataSplitTable[,i]), sp.names.ok])
      speciesData.evaluation[,1:sum(!DataSplitTable[,i]),i] <- t(sp.data[which(!DataSplitTable[,i]), sp.names.ok])
    }
    
    #Creating the arrays to save the results
    singleSpecies.ensembleEvaluationScore <- array(data=NA, dim=c(length(eval.metrics.names), length(sp.names.ok), NbRunEval), dimnames = list(eval.metrics.names,sp.names.ok,1:NbRunEval))
    singleSpecies.ensembleVariableImportance <- array(data=NA, dim=c(NbPredictors, length(sp.names.ok), NbRunEval), dimnames = list(NamesPredictors,sp.names.ok,1:NbRunEval))
    singleSpecies.calibrationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(DataSplitTable)),  NbRunEval), dimnames=list(sp.names.ok, paste("c",1:max(colSums(DataSplitTable)),sep="_"), 1:NbRunEval))
    singleSpecies.evaluationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(!DataSplitTable)), NbRunEval), dimnames=list(sp.names.ok, paste("e",1:max(colSums(!DataSplitTable)),sep="_"), 1:NbRunEval))
    
    #Running the models for all species and cross-validation runs| This part needs to be extracted if run on a cluster to increase efficiency
    if(parallel){
      sfInit(parallel=TRUE, cpus=cpus)
      sfLibrary('biomod2', character.only=TRUE)
      sfLapply(sp.names.ok, 
               BiomodSF,
               DataSplitTable=DataSplitTable, 
               sp.data=sp.data, 
               env.data=env.data, 
               xy=xy, 
               models=models.sdm, 
               models.options=modeling.options.sdm, 
               eval.metrics=eval.metrics.sdm, 
               ensemble.metric=ensemble.metric.sdm,
               VarImport = VarImport)
      sfStop( nostop=FALSE )
    }else{
      lapply(sp.names.ok, 
             BiomodSF,
             DataSplitTable=DataSplitTable, 
             sp.data=sp.data, 
             env.data=env.data, 
             xy=xy, 
             models=models.sdm, 
             models.options=modeling.options.sdm, 
             eval.metrics=eval.metrics.sdm, 
             ensemble.metric=ensemble.metric.sdm,
             VarImport = VarImport)
    }
    
    #Gathering the results
    
    #For the SDMs
    for(i in sp.names.ok){
      load(paste(i,"/",i,".ccvensemble.models.out", sep=""))
      
      #Single model evaluation
      temp.evaluations <- get_evaluations(eval(parse(text=paste(i,".ccvensemble.models.out",sep=""))))
      for(l in 1:length(temp.evaluations)){
        singleSpecies.ensembleEvaluationScore[,i,l] <- temp.evaluations[[l]][,1]
      }
      
      #Single model variable importance
      temp.variableimprtance <- get_variables_importance(eval(parse(text=paste(i,".ccvensemble.models.out",sep=""))))
      singleSpecies.ensembleVariableImportance[,i,] <- round(apply(temp.variableimprtance,c(1,3), mean, na.rm = TRUE),2)
      
      #Single model predictions
      temp.predictions <- get_predictions(eval(parse(text=paste(i,".ccvensemble.models.out",sep=""))))
      for(l in 1:dim(temp.predictions)[2]){
        singleSpecies.calibrationSites.ensemblePredictions[i,1:sum(DataSplitTable[,l]),l] <- temp.predictions[which(DataSplitTable[,l]),l]
        singleSpecies.evaluationSites.ensemblePredictions[i,1:sum(!DataSplitTable[,l]),l] <- temp.predictions[which(!DataSplitTable[,l]),l]
      }
    }
    
  }
  
  #Running the stuff with ONLY ESMs #################################################################################################
  if(ESM=="ALL"){
    
    #Selecting the species than can be run (enough data)
    sp.names.all <- rownames(SpRunMatrix)
    sp.names.ok <- intersect(names(which(apply(SpRunMatrix,1,min) >= minNbPredictors*2)), names(which(apply(min(colSums(DataSplitTable))-SpRunMatrix,1,min) >= minNbPredictors*2)))
    sp.names.droped <- setdiff(sp.names.all, sp.names.ok) 
    message(paste("The following species will not be modelled due to limited presence data: ", paste(sp.names.droped, sep="", collapse=", "),"\n\n",sep=""))
    
    #Creating the speciesData.calibration and speciesData.evaluation
    speciesData.calibration <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(DataSplitTable)),  NbRunEval), dimnames=list(sp.names.ok, paste("c",1:max(colSums(DataSplitTable)),sep="_"), 1:NbRunEval))
    speciesData.evaluation <- singleSpecies.evaluationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(!DataSplitTable)), NbRunEval), dimnames=list(sp.names.ok, paste("e",1:max(colSums(!DataSplitTable)),sep="_"), 1:NbRunEval))
    for(i in 1:NbRunEval){
      speciesData.calibration[,1:sum(DataSplitTable[,i]),i] <- t(sp.data[which(DataSplitTable[,i]), sp.names.ok])
      speciesData.evaluation[,1:sum(!DataSplitTable[,i]),i] <- t(sp.data[which(!DataSplitTable[,i]), sp.names.ok])
    }
    
    #Creating the arrays to save the results
    singleSpecies.ensembleEvaluationScore <- array(data=NA, dim=c(length(eval.metrics.names), length(sp.names.ok), NbRunEval), dimnames = list(eval.metrics.names,sp.names.ok,1:NbRunEval))
    singleSpecies.ensembleVariableImportance <- array(data=NA, dim=c(NbPredictors, length(sp.names.ok), NbRunEval), dimnames = list(NamesPredictors,sp.names.ok,1:NbRunEval))
    singleSpecies.calibrationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(DataSplitTable)),  NbRunEval), dimnames=list(sp.names.ok, paste("c",1:max(colSums(DataSplitTable)),sep="_"), 1:NbRunEval))
    singleSpecies.evaluationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(!DataSplitTable)), NbRunEval), dimnames=list(sp.names.ok, paste("e",1:max(colSums(!DataSplitTable)),sep="_"), 1:NbRunEval))
    
    #Running the models for all species and cross-validation runs| This part needs to be extracted if run on a cluster to increase efficiency
    if(parallel){
      sfInit(parallel=TRUE, cpus=cpus)
      sfLibrary('biomod2', character.only=TRUE)
      sfLibrary('ecospat', character.only=TRUE)
      sfLibrary('gtools', character.only=TRUE)
      sfLapply(sp.names.ok, 
               ESMSF, 
               DataSplitTable=DataSplitTable, 
               sp.data=sp.data, 
               env.data=env.data, 
               xy=xy, 
               models=models.esm, 
               models.options=modeling.options.esm, 
               ensemble.metric=ensemble.metric.esm)
      sfStop( nostop=FALSE)
    }else{
      lapply(sp.names.ok,
             ESMSF, 
             DataSplitTable=DataSplitTable, 
             sp.data=sp.data, 
             env.data=env.data, 
             xy=xy, 
             models=models.esm, 
             models.options=modeling.options.esm, 
             ensemble.metric=ensemble.metric.esm)
    }
    
    #Gathering the results
    #For the ESMs
    for(i in sp.names.ok){
      load(list.files(path=paste("ESM.BIOMOD.output_",i,sep=""), pattern="ESM_EnsembleModeling", full.names = TRUE))
      output_EF <- eval(parse(text="output"))
      load(list.files(path=paste("ESM.BIOMOD.output_",i,sep=""), pattern="ESM_Modeling", full.names = TRUE))
      
      #Single model evaluation
      singleSpecies.ensembleEvaluationScore[,i,] <- t(output_EF$ESM.evaluations[seq(ef.counter,dim(output_EF$ESM.evaluations)[1], ef.counter),c(5,11,6)])
      
      #Single model variable importance
      singleSpecies.ensembleVariableImportance[,i,] <- round(get.ESMvariableContribution(output_EF = output_EF, output = eval(parse(text="output")), NamesPredictors = NamesPredictors),2)
      
      #Single model predictions
      for(l in 1:NbRunEval){
        singleSpecies.calibrationSites.ensemblePredictions[i,1:sum(DataSplitTable[,l]),l] <- output_EF$ESM.fit[which(DataSplitTable[,l]),ef.counter*l+1]
        singleSpecies.evaluationSites.ensemblePredictions[i,1:sum(!DataSplitTable[,l]),l] <- output_EF$ESM.fit[which(!DataSplitTable[,l]),ef.counter*l+1]
      }
    }
  }
  #Running the with biomod2 and ESMs #################################################################################################
  if(ESM=="YES"){
    
    #Selecting the species than can be run (enough data)
    sp.names.all <- rownames(SpRunMatrix)
    sp.names.bm.ok <- intersect(names(which(apply(SpRunMatrix,1, min) >= minNbPredictors*NbPredictors)), names(which(apply(min(colSums(DataSplitTable))-SpRunMatrix,1,min) >= minNbPredictors*NbPredictors)))
    message(paste("The following species will be run with standard biomod2 models: ", paste(sp.names.bm.ok, sep="", collapse=", "),"\n\n",sep=""))
    sp.names.bm.droped <- setdiff(sp.names.all, sp.names.bm.ok)
    if(length(sp.names.bm.droped)>1){
      sp.names.esm.ok <- intersect(names(which(apply(SpRunMatrix[sp.names.bm.droped,],1,min) >= minNbPredictors*2)), names(which(apply(min(colSums(DataSplitTable))-SpRunMatrix[sp.names.bm.droped,],1,min) >= minNbPredictors*2)))
      message(paste("The following species will be run with ESM models: ", paste(sp.names.esm.ok, sep="", collapse=", "),"\n\n",sep=""))
      sp.names.droped <- setdiff(sp.names.all, c(sp.names.bm.ok, sp.names.esm.ok))
      message(paste("The following species will not be modelled due to limited presence data: ", paste(sp.names.droped, sep="", collapse=", "),"\n\n",sep=""))
      sp.names.ok <- sort(c(sp.names.bm.ok, sp.names.esm.ok))
    }else{
      if(length(sp.names.bm.droped==1)){
        if(min(SpRunMatrix[sp.names.bm.droped,]) >= minNbPredictors*2 & min(colSums(DataSplitTable))-min(SpRunMatrix[sp.names.bm.droped,])>= minNbPredictors*2)
        sp.names.esm.ok <- sp.names.bm.droped
        message(paste("The following species will be run with ESM models: ", paste(sp.names.esm.ok, sep="", collapse=", "),"\n\n",sep=""))
        message(paste("The following species will not be modelled due to limited presence data:","\n\n", sep=""))
        sp.names.ok <- sort(c(sp.names.bm.ok, sp.names.esm.ok))
      }else{
        sp.names.esm.ok <- NULL
        message(paste("The following species will be run with ESM models:","\n\n", sep=""))
        message(paste("The following species will not be modelled due to limited presence data:","\n\n", sep=""))
        sp.names.ok <- sort(sp.names.bm.ok)
      }
    }
    
    
    #Creating the speciesData.calibration and speciesData.evaluation
    speciesData.calibration <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(DataSplitTable)),  NbRunEval), dimnames=list(sp.names.ok, paste("c",1:max(colSums(DataSplitTable)),sep="_"), 1:NbRunEval))
    speciesData.evaluation <- singleSpecies.evaluationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(!DataSplitTable)), NbRunEval), dimnames=list(sp.names.ok, paste("e",1:max(colSums(!DataSplitTable)),sep="_"), 1:NbRunEval))
    for(i in 1:NbRunEval){
      speciesData.calibration[,1:sum(DataSplitTable[,i]),i] <- t(sp.data[which(DataSplitTable[,i]), sp.names.ok])
      speciesData.evaluation[,1:sum(!DataSplitTable[,i]),i] <- t(sp.data[which(!DataSplitTable[,i]), sp.names.ok])
    }
    
    #Creating the arrays to save the results
    singleSpecies.ensembleEvaluationScore <- array(data=NA, dim=c(length(eval.metrics.names), length(sp.names.ok), NbRunEval), dimnames = list(eval.metrics.names,sp.names.ok,1:NbRunEval))
    singleSpecies.ensembleVariableImportance <- array(data=NA, dim=c(NbPredictors, length(sp.names.ok), NbRunEval), dimnames = list(NamesPredictors,sp.names.ok,1:NbRunEval))
    singleSpecies.calibrationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(DataSplitTable)),  NbRunEval), dimnames=list(sp.names.ok, paste("c",1:max(colSums(DataSplitTable)),sep="_"), 1:NbRunEval))
    singleSpecies.evaluationSites.ensemblePredictions <- array(data=NA, dim=c(length(sp.names.ok), max(colSums(!DataSplitTable)), NbRunEval), dimnames=list(sp.names.ok, paste("e",1:max(colSums(!DataSplitTable)),sep="_"), 1:NbRunEval))
    
    #Running the models for all species and cross-validation runs| This part needs to be extracted if run on a cluster to increase efficiency
    if(parallel){
      sfInit(parallel=TRUE, cpus=cpus)
      sfLibrary('biomod2', character.only=TRUE)
      sfLibrary('ecospat', character.only=TRUE)
      sfLibrary('gtools', character.only=TRUE)
      sfLapply(sp.names.bm.ok, 
               BiomodSF,
               DataSplitTable=DataSplitTable, 
               sp.data=sp.data, 
               env.data=env.data, 
               xy=xy, models=models.sdm, 
               models.options=modeling.options.sdm, 
               eval.metrics=eval.metrics.sdm, 
               ensemble.metric=ensemble.metric.sdm,
               VarImport = VarImport)
      sfLapply(sp.names.esm.ok, 
               ESMSF, 
               DataSplitTable=DataSplitTable, 
               sp.data=sp.data, 
               env.data=env.data, 
               xy=xy, 
               models=models.esm, 
               models.options=modeling.options.esm, 
               ensemble.metric=ensemble.metric.esm)
      sfStop( nostop=FALSE )
    }else{
      lapply(sp.names.bm.ok, 
             BiomodSF,
             DataSplitTable=DataSplitTable, 
             sp.data=sp.data, 
             env.data=env.data, 
             xy=xy, models=models.sdm, 
             models.options=modeling.options.sdm, 
             eval.metrics=eval.metrics.sdm, 
             ensemble.metric=ensemble.metric.sdm,
             VarImport = VarImport)
      lapply(sp.names.esm.ok, 
             ESMSF, 
             DataSplitTable=DataSplitTable, 
             sp.data=sp.data, 
             env.data=env.data, 
             xy=xy, 
             models=models.esm, 
             models.options=modeling.options.esm, 
             ensemble.metric=ensemble.metric.esm)
    }
    
    #Gathering the results
    
    #For the SDMs
    for(i in sp.names.bm.ok){
      load(paste(i,"/",i,".ccvensemble.models.out", sep=""))
      
      #Single model evaluation
      temp.evaluations <- get_evaluations(eval(parse(text=paste(i,".ccvensemble.models.out",sep=""))))
      for(l in 1:length(temp.evaluations)){
        singleSpecies.ensembleEvaluationScore[,i,l] <- temp.evaluations[[l]][,1]
      }
      
      #Single model variable importance
      temp.variableimprtance <- get_variables_importance(eval(parse(text=paste(i,".ccvensemble.models.out",sep=""))))
      singleSpecies.ensembleVariableImportance[,i,] <- round(apply(temp.variableimprtance,c(1,3), mean, na.rm = TRUE),2)
      
      #Single model predictions
      temp.predictions <- get_predictions(eval(parse(text=paste(i,".ccvensemble.models.out",sep=""))))
      for(l in 1:dim(temp.predictions)[2]){
        singleSpecies.calibrationSites.ensemblePredictions[i,1:sum(DataSplitTable[,l]),l] <- temp.predictions[which(DataSplitTable[,l]),l]
        singleSpecies.evaluationSites.ensemblePredictions[i,1:sum(!DataSplitTable[,l]),l] <- temp.predictions[which(!DataSplitTable[,l]),l]
      }
    }
    
    #For the ESMs
    if(!is.null(sp.names.esm.ok)){
      for(i in sp.names.esm.ok){
        load(list.files(path=paste("ESM.BIOMOD.output_",i,sep=""), pattern="ESM_EnsembleModeling", full.names = TRUE))
        output_EF <- eval(parse(text="output"))
        load(list.files(path=paste("ESM.BIOMOD.output_",i,sep=""), pattern="ESM_Modeling", full.names = TRUE))
        
        #Single model evaluation
        singleSpecies.ensembleEvaluationScore[,i,] <- t(output_EF$ESM.evaluations[seq(ef.counter,dim(output_EF$ESM.evaluations)[1], ef.counter),c(5,11,6)])
        
        #Single model variable importance
        singleSpecies.ensembleVariableImportance[,i,] <- round(get.ESMvariableContribution(output_EF = output_EF, output = eval(parse(text="output")), NamesPredictors = NamesPredictors),2)
        
        #Single model predictions
        for(l in 1:NbRunEval){
          singleSpecies.calibrationSites.ensemblePredictions[i,1:sum(DataSplitTable[,l]),l] <- output_EF$ESM.fit[which(DataSplitTable[,l]),ef.counter*l+1]
          singleSpecies.evaluationSites.ensemblePredictions[i,1:sum(!DataSplitTable[,l]),l] <- output_EF$ESM.fit[which(!DataSplitTable[,l]),ef.counter*l+1]
        }
      }
    }
  }
  
  #Creating the average probability predictions for each site  ############################################################################
  all.predictions.caliSites <- array(data=NA, 
                                     dim=c(dim(sp.data),dim(DataSplitTable)[2]),
                                     dimnames=list(unlist(dimnames(sp.data)[1]),
                                                   unlist(dimnames(sp.data)[2]),
                                                   1:dim(DataSplitTable)[2]))
  all.predictions.evalSites <- array(data=NA, 
                                     dim=c(dim(sp.data),dim(DataSplitTable)[2]),
                                     dimnames=list(unlist(dimnames(sp.data)[1]),
                                                   unlist(dimnames(sp.data)[2]),
                                                   1:dim(DataSplitTable)[2]))
  
  #Extracting the data from the CCV.Modelling output
  for(i in 1:dim(DataSplitTable)[2]){
    all.predictions.caliSites[DataSplitTable[,i],,i] <- t(singleSpecies.calibrationSites.ensemblePredictions[,,i])[1:dim(all.predictions.caliSites[DataSplitTable[,i],,i])[1]]
    all.predictions.evalSites[!DataSplitTable[,i],,i] <- t(singleSpecies.evaluationSites.ensemblePredictions[,,i])[1:dim(all.predictions.evalSites[!DataSplitTable[,i],,i])[1]]
  }
  
  #Making the average prediction per site
  allSites.averagePredictions.cali <- apply(all.predictions.caliSites, 1:2, mean, na.rm = TRUE)
  allSites.averagePredictions.eval <- apply(all.predictions.evalSites, 1:2, mean, na.rm = TRUE)
  
  #Writing the final output files ##############################################################################################
  save(singleSpecies.ensembleEvaluationScore, file="singleSpecies.ensembleEvaluationScore.RData")
  save(singleSpecies.calibrationSites.ensemblePredictions, file="singleSpecies.calibrationSites.ensemblePredictions.RData")
  save(singleSpecies.evaluationSites.ensemblePredictions, file="singleSpecies.evaluationSites.ensemblePredictions.RData")
  save(singleSpecies.ensembleVariableImportance, file="singleSpecies.ensembleVariableImportance.RData")
  save(DataSplitTable, file="DataSplitTable.RData")
  save(speciesData.calibration, file="speciesData.calibration.RData")
  save(speciesData.evaluation, file="speciesData.evaluation.RData")
  save(allSites.averagePredictions.cali, file="allSites.averagePredictions.cali.RData")
  save(allSites.averagePredictions.eval, file="allSites.averagePredictions.eval.RData")
  ccv.modeling.data <- list(modeling.id = modeling.id,
                            output.files = c("singleSpecies.ensembleEvaluationScore.RData",
                                             "singleSpecies.calibrationSites.ensemblePredictions.RData",
                                             "singleSpecies.evaluationSites.ensemblePredictions.RData",
                                             "singleSpecies.ensembleVariableImportance.RData",
                                             "DataSplitTable.RData",
                                             "speciesData.calibration.RData",
                                             "speciesData.evaluation.RData",
                                             "allSites.averagePredictions.cali.RData",
                                             "allSites.averagePredictions.eval.RData"),
                            speciesData.calibration = speciesData.calibration,
                            speciesData.evaluation = speciesData.evaluation,
                            speciesData.full = sp.data,
                            DataSplitTable = DataSplitTable,
                            singleSpecies.ensembleEvaluationScore = singleSpecies.ensembleEvaluationScore,
                            singleSpecies.ensembleVariableImportance = singleSpecies.ensembleVariableImportance,
                            singleSpecies.calibrationSites.ensemblePredictions=singleSpecies.calibrationSites.ensemblePredictions,
                            singleSpecies.evaluationSites.ensemblePredictions=singleSpecies.evaluationSites.ensemblePredictions,
                            allSites.averagePredictions.cali=allSites.averagePredictions.cali,
                            allSites.averagePredictions.eval=allSites.averagePredictions.eval)
  save(ccv.modeling.data, file=paste("../",modeling.id,".ccv.modeling.RData", sep=""))
  setwd("../")
  return(ccv.modeling.data)
}


#####################################################
# c) ecospat.CCV.communityEvaluation.bin
#####################################################

#DESCRIPTION
# Evaluates the community predictions based on different thresholding techniques. 
# This function binarised the individual predictions and then stacks the binary maps. There are a range
# of different binarisation methods available and the resulting species assemblages get evaluated 
# by different community evaluation metrics


#FUNCTION'S ARGUMENTS
#ccv.modeling.data     an output from ecospat.CCV.modeling function
#thresholds            a selection of thresholding methods used for the community building (FIXED, MAX.KAPPA, MAX.ACCURACY, MAX.TSS, SENS_SPEC, MAX.ROC, OBS.PREVALENCE, AVG.PROBABILITY, MCE, PS_SDM, MEM)
#community.metrics     a selection of community evaluation metrics to be calculated for each selected threshold (SR.deviation, community.AUC, community.overprediction ,community.underprediction, community.accuracy, community.sensitivity, community.specificity, community.kappa, community.tss, Sorensen, Jaccard, Simpson)
#parallel              binary variable if parallel computing is allowed
#cpus                  if parallel true the number of cpus to use
#fix.threshold         if FIXED selected as threshold this is the value used
#MCE                   if MCE is selected as threshold this is the maximal commison error
#MEM                   if MEM is selected as threshold these are the SR predictions used


#DETAILS


#VALUES
#DataSplitTable                                           a matrix with TRUE/FALSE for each model run (TRUE=Calibration point, FALSE=Evaluation point)
#CommunityEvaluationMetrics.CalibrationSites              a 4-dimensional array containing the community evaluation metrics for the calibartion sites of each run (NA means that the site was used for evaluation)
#CommunityEvaluationMetrics.EvaluationSites               a 4-dimensional array containing the community evaluation metrics for the evaluation sites of each run (NA means that the site was used for calibaration)
#PA.allSites                                              a 4-dimensional array of the binary prediction for all sites and runs under the different thresholding appraoches.



#AUTHORS
# Daniel Scherrer <daniel.j.a.scherrer@gmail.com>


#REFERENCES
# Scherrer, D., D'Amen, M., Mateo, M.R.G., Fernandes, R.F. & Guisan , A. (2018) 
#   How to best threshold and validate stacked species assemblages? Community optimisation might hold the answer. 
#   Methods in Ecology and Evolution, 9(10): 2155-2166.


#SEE ALSO
# ecospat.CCV.modelling


ecospat.CCV.communityEvaluation.bin <- function(ccv.modeling.data,
                                                thresholds= c("MAX.KAPPA", "MAX.ROC","PS_SDM"),
                                                community.metrics=c("SR.deviation","Sorensen"),
                                                parallel=FALSE,
                                                cpus=4,
                                                fix.threshold=0.5,
                                                MCE=5,
                                                MEM=NULL){

    #Loading the packages needed (they should all be installed by ecospat library)
  #require(snowfall)
  #require(PresenceAbsence)
  
  #Checking all the input data
  stopifnot(names(ccv.modeling.data)==c("modeling.id",
                                        "output.files",
                                        "speciesData.calibration",
                                        "speciesData.evaluation",
                                        "speciesData.full",
                                        "DataSplitTable",
                                        "singleSpecies.ensembleEvaluationScore",
                                        "singleSpecies.ensembleVariableImportance",
                                        "singleSpecies.calibrationSites.ensemblePredictions",
                                        "singleSpecies.evaluationSites.ensemblePredictions",
                                        "allSites.averagePredictions.cali",
                                        "allSites.averagePredictions.eval"))
  possible.thresholds <- c("FIXED", 
                           "MAX.KAPPA",
                           "MAX.ACCURACY",
                           "MAX.TSS",
                           "SENS_SPEC",
                           "MAX.ROC",
                           "OBS.PREVALENCE",
                           "AVG.PROBABILITY",
                           "MCE",
                           "PS_SDM",
                           "MEM")
  stopifnot(thresholds %in% possible.thresholds)
  stopifnot(community.metrics %in% c("SR.deviation",
                                     "community.overprediction",
                                     "community.underprediction",
                                     "community.accuracy",
                                     "community.sensitivity",
                                     "community.specificity",
                                     "community.kappa",
                                     "community.tss",
                                     "Sorensen",
                                     "Jaccard",
                                     "Simpson"))
  stopifnot(is.logical(parallel))
  stopifnot(cpus>=1)
  stopifnot(!("FIXED" %in% thresholds & (fix.threshold<=0 | fix.threshold>=1)))
  stopifnot(!("MCE" %in% thresholds & (MCE<=0 | MCE>=100)))
  stopifnot(!("MEM" %in% thresholds & length(MEM)!=dim(ccv.modeling.data$speciesData.full)[1])) 
  
  
  #############################################################################################################
  #SUBFUNCTIONS################################################################################################
  #############################################################################################################
  
  #Calculating the community evaluation metrics
  community.metrics.calculation <- function(errors, potential.community.metrics){
    temp.matrix <- matrix(data=NA, nrow=1, ncol=length(potential.community.metrics))
    a <- length(which(errors == 3)) #True Positives
    b <- length(which(errors == 2)) #False Positives
    c <- length(which(errors == 1)) #False Negatives
    d <- length(which(errors == 0)) #True Negatives
    n <- a+b+c+d
    
    #SR.deviance
    temp.matrix[1] <- b-c
    
    #Community.overprediction
    if(b==0 & d==0){
      temp.matrix[2] <- 0
    }else{
      temp.matrix[2] <-  round(b/(b + d), digits=3)
    }
    
    #Community.underprediction
    if(a==0 & c==0){
      temp.matrix[3] <- 0
    }else{
      temp.matrix[3] <- round(c/(a + c), digits=3)
    }
    
    #Community.accuracy
    if(n==0){
      temp.matrix[4] <- 1
    }else{
      temp.matrix[4] <- round((a + d)/n, digits=3)
    }
    
    #Community.sensitivity
    if(a==0 & c==0){
      temp.matrix[5] <- 1
    }else{
      temp.matrix[5] <- round(a/(a + c), digits=3)
    }
    
    #Community.specificity
    if(b==0 & d==0){
      temp.matrix[6] <- 1
    }else{
      temp.matrix[6] <- round(d/(b + d), digits=3)
    }
    
    #community.kappa
    if(n==0){
      temp.matrix[7] <- 1
    }else{
      temp.matrix[7] <- round((((a + d)/n) - (((a + c) * (a + b) + (b + d) * (d + c))/(n^2)))/(1 - (((a + c) * (a + b) + (b + d) * (d + c))/(n^2))), digits=3)
    }
    
    #community.tss
    temp.matrix[8] <- round(temp.matrix[5] + temp.matrix[6] - 1, digits=3)
    
    #Sorensen
    if(a==0 & b==0 & c==0){
      temp.matrix[9] <- 1
    }else{
      temp.matrix[9] <- round((2 * a)/(2 * a + b + c), digits=3)
    }
    
    #Jaccard
    if(a==0 & b==0 & c==0){
      temp.matrix[10] <- 1
    }else{
      temp.matrix[10] <- round(a/(a+b+c), digits=3)
    }
    
    #Simpson
    if((a==0 & b==0) | (a==0 & c==0)){
      temp.matrix[11] <- 1
    }else{
      temp.matrix[11] <- round(a/min(c(a+b,a+c)), digits=3)
    }
    
    return(temp.matrix)
  }
  
  #Community evaluation metrics calculation
  community.compairison <- function(sp.data.cali, sp.data.eval, PA.cali, PA.eval, community.metrics, save.dir, run){
    
    #Naming the potential community evalutaion metrics
    potential.community.metrics <- c("SR.deviation",
                                     "community.overprediction",
                                     "community.underprediction",
                                     "community.accuracy",
                                     "community.sensitivity",
                                     "community.specificity",
                                     "community.kappa",
                                     "community.tss",
                                     "Sorensen",
                                     "Jaccard",
                                     "Simpson")
    #Arrays to save the metrics
    community.metrics.cali <- array(data=NA, 
                                    dim=c(dim(sp.data.cali)[1],length(potential.community.metrics)), 
                                    dimnames=list(unlist(dimnames(sp.data.cali)[1]),potential.community.metrics))
    community.metrics.eval <- array(data=NA, 
                                    dim=c(dim(sp.data.eval)[1],length(potential.community.metrics)), 
                                    dimnames=list(unlist(dimnames(sp.data.eval)[1]),potential.community.metrics))
    
    #Calculation of the error matrix [Congencenty table]
    error.matrix.cali <- 2 * PA.cali + sp.data.cali
    error.matrix.eval <- 2 * PA.eval + sp.data.eval
    
    #Calculating all the metrics and then selection of the ones desired for the output
    community.metrics.cali[,] <- t(apply(error.matrix.cali, 1, community.metrics.calculation, potential.community.metrics=potential.community.metrics))
    community.metrics.cali.selected <- community.metrics.cali[,community.metrics]
    community.metrics.eval[,] <- t(apply(error.matrix.eval, 1, community.metrics.calculation, potential.community.metrics=potential.community.metrics))
    community.metrics.eval.selected <- community.metrics.eval[,community.metrics]
    
    #Writing all the results
    save(community.metrics.cali.selected, file=paste(save.dir,"community.metrics.cali_",run,".RData", sep=""))
    save(community.metrics.eval.selected, file=paste(save.dir,"community.metrics.eval_",run,".RData", sep=""))
  }
  
  
  #Calculating the community evaluations ######################################################################
  community.thresholding <- function(run, ccv.modeling.data, thresholds, community.metrics, fix.threshold, MCE, MEM){
    
    #-------------------------------------------------
    #Thresholding by fixed threshold
    #-------------------------------------------------
    if("FIXED" %in% thresholds){
      #Converting probabilities to binary
      #Calibration data
      PA.FIXED.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.FIXED.cali[PA.FIXED.cali >= fix.threshold] <- 1
      PA.FIXED.cali[PA.FIXED.cali < fix.threshold] <- 0
      #Evaluation data
      PA.FIXED.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      PA.FIXED.eval[PA.FIXED.eval >= fix.threshold] <- 1
      PA.FIXED.eval[PA.FIXED.eval < fix.threshold] <- 0
      
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/FIXED", sep=""), recursive=TRUE)
      save(PA.FIXED.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/FIXED/PA.FIXED.cali_",run,".RData", sep=""))
      save(PA.FIXED.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/FIXED/PA.FIXED.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.FIXED.cali,
                            PA.eval = PA.FIXED.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/FIXED/FIXED_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by Max Kappa
    #-------------------------------------------------
    if("MAX.KAPPA" %in% thresholds){
      PA.MAX.KAPPA.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.MAX.KAPPA.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.MAX.KAPPA.cali)[2]){
        if(sum(!is.na(PA.MAX.KAPPA.cali[,s]))==0){
          PA.MAX.KAPPA.cali[,s] <- NA
          PA.MAX.KAPPA.eval[,s] <- NA
        }else{
          MAX.KAPPA.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.MAX.KAPPA.cali)[1]),
                                                                            t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                            t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                    threshold=101, opt.methods=4)[,2]
          PA.MAX.KAPPA.cali[PA.MAX.KAPPA.cali[,s] >= MAX.KAPPA.threshold,s] <- 1
          PA.MAX.KAPPA.cali[PA.MAX.KAPPA.cali[,s] <= MAX.KAPPA.threshold,s] <- 0
          PA.MAX.KAPPA.eval[PA.MAX.KAPPA.eval[,s] >= MAX.KAPPA.threshold,s] <- 1
          PA.MAX.KAPPA.eval[PA.MAX.KAPPA.eval[,s] <= MAX.KAPPA.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.KAPPA", sep=""), recursive=TRUE)
      save(PA.MAX.KAPPA.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.KAPPA/PA.MAX.KAPPA.cali_",run,".RData", sep=""))
      save(PA.MAX.KAPPA.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.KAPPA/PA.MAX.KAPPA.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.MAX.KAPPA.cali,
                            PA.eval = PA.MAX.KAPPA.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.KAPPA/MAX.KAPPA_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by Max Accuracy
    #-------------------------------------------------
    if("MAX.ACCURACY" %in% thresholds){
      PA.MAX.ACCURACY.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.MAX.ACCURACY.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.MAX.ACCURACY.cali)[2]){
        if(sum(!is.na(PA.MAX.ACCURACY.cali[,s]))==0){
          PA.MAX.ACCURACY.cali[,s] <- NA
          PA.MAX.ACCURACY.eval[,s] <- NA
        }else{
          MAX.ACCURACY.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.MAX.ACCURACY.cali)[1]),
                                                                               t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                               t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                       threshold=101, opt.methods=5)[,2]
          PA.MAX.ACCURACY.cali[PA.MAX.ACCURACY.cali[,s] >= MAX.ACCURACY.threshold,s] <- 1
          PA.MAX.ACCURACY.cali[PA.MAX.ACCURACY.cali[,s] <= MAX.ACCURACY.threshold,s] <- 0
          PA.MAX.ACCURACY.eval[PA.MAX.ACCURACY.eval[,s] >= MAX.ACCURACY.threshold,s] <- 1
          PA.MAX.ACCURACY.eval[PA.MAX.ACCURACY.eval[,s] <= MAX.ACCURACY.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ACCURACY", sep=""), recursive=TRUE)
      save(PA.MAX.ACCURACY.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ACCURACY/PA.MAX.ACCURACY.cali_",run,".RData", sep=""))
      save(PA.MAX.ACCURACY.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ACCURACY/PA.MAX.ACCURACY.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.MAX.ACCURACY.cali,
                            PA.eval = PA.MAX.ACCURACY.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ACCURACY/MAX.ACCURACY_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by Max TSS
    #-------------------------------------------------
    if("MAX.TSS" %in% thresholds){
      PA.MAX.TSS.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.MAX.TSS.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.MAX.TSS.cali)[2]){
        if(sum(!is.na(PA.MAX.TSS.cali[,s]))==0){
          PA.MAX.TSS.cali[,s] <- NA
          PA.MAX.TSS.eval[,s] <- NA
        }else{
          MAX.TSS.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.MAX.TSS.cali)[1]),
                                                                          t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                          t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                  threshold=101, opt.methods=3)[,2]
          PA.MAX.TSS.cali[PA.MAX.TSS.cali[,s] >= MAX.TSS.threshold,s] <- 1
          PA.MAX.TSS.cali[PA.MAX.TSS.cali[,s] <= MAX.TSS.threshold,s] <- 0
          PA.MAX.TSS.eval[PA.MAX.TSS.eval[,s] >= MAX.TSS.threshold,s] <- 1
          PA.MAX.TSS.eval[PA.MAX.TSS.eval[,s] <= MAX.TSS.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.TSS", sep=""), recursive=TRUE)
      save(PA.MAX.TSS.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.TSS/PA.MAX.TSS.cali_",run,".RData", sep=""))
      save(PA.MAX.TSS.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.TSS/PA.MAX.TSS.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.MAX.TSS.cali,
                            PA.eval = PA.MAX.TSS.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.TSS/MAX.TSS_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by SENS_SPEC
    #-------------------------------------------------
    if("SENS_SPEC" %in% thresholds){
      PA.SENS_SPEC.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.SENS_SPEC.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.SENS_SPEC.cali)[2]){
        if(sum(!is.na(PA.SENS_SPEC.cali[,s]))==0){
          PA.SENS_SPEC.cali[,s] <- NA
          PA.SENS_SPEC.eval[,s] <- NA
        }else{
          SENS_SPEC.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.SENS_SPEC.cali)[1]),
                                                                            t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                            t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                    threshold=101, opt.methods=2)[,2]
          PA.SENS_SPEC.cali[PA.SENS_SPEC.cali[,s] >= SENS_SPEC.threshold,s] <- 1
          PA.SENS_SPEC.cali[PA.SENS_SPEC.cali[,s] <= SENS_SPEC.threshold,s] <- 0
          PA.SENS_SPEC.eval[PA.SENS_SPEC.eval[,s] >= SENS_SPEC.threshold,s] <- 1
          PA.SENS_SPEC.eval[PA.SENS_SPEC.eval[,s] <= SENS_SPEC.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/SENS_SPEC", sep=""), recursive=TRUE)
      save(PA.SENS_SPEC.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/SENS_SPEC/PA.SENS_SPEC.cali_",run,".RData", sep=""))
      save(PA.SENS_SPEC.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/SENS_SPEC/PA.SENS_SPEC.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.SENS_SPEC.cali,
                            PA.eval = PA.SENS_SPEC.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/SENS_SPEC/SENS_SPEC_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by Max ROC
    #-------------------------------------------------
    if("MAX.ROC" %in% thresholds){
      PA.MAX.ROC.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.MAX.ROC.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.MAX.ROC.cali)[2]){
        if(sum(!is.na(PA.MAX.ROC.cali[,s]))==0){
          PA.MAX.ROC.cali[,s] <- NA
          PA.MAX.ROC.eval[,s] <- NA
        }else{
          MAX.ROC.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.MAX.ROC.cali)[1]),
                                                                          t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                          t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                  threshold=101, opt.methods=9)[,2]
          PA.MAX.ROC.cali[PA.MAX.ROC.cali[,s] >= MAX.ROC.threshold,s] <- 1
          PA.MAX.ROC.cali[PA.MAX.ROC.cali[,s] <= MAX.ROC.threshold,s] <- 0
          PA.MAX.ROC.eval[PA.MAX.ROC.eval[,s] >= MAX.ROC.threshold,s] <- 1
          PA.MAX.ROC.eval[PA.MAX.ROC.eval[,s] <= MAX.ROC.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ROC", sep=""), recursive=TRUE)
      save(PA.MAX.ROC.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ROC/PA.MAX.ROC.cali_",run,".RData", sep=""))
      save(PA.MAX.ROC.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ROC/PA.MAX.ROC.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.MAX.ROC.cali,
                            PA.eval = PA.MAX.ROC.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/MAX.ROC/MAX.ROC_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by Observed Prevalence
    #-------------------------------------------------
    if("OBS.PREVALENCE" %in% thresholds){
      PA.OBS.PREVALENCE.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.OBS.PREVALENCE.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.OBS.PREVALENCE.cali)[2]){
        if(sum(!is.na(PA.OBS.PREVALENCE.cali[,s]))==0){
          PA.OBS.PREVALENCE.cali[,s] <- NA
          PA.OBS.PREVALENCE.eval[,s] <- NA
        }else{
          OBS.PREVALENCE.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.OBS.PREVALENCE.cali)[1]),
                                                                                 t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                                 t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                         threshold=101, opt.methods=6)[,2]
          PA.OBS.PREVALENCE.cali[PA.OBS.PREVALENCE.cali[,s] >= OBS.PREVALENCE.threshold,s] <- 1
          PA.OBS.PREVALENCE.cali[PA.OBS.PREVALENCE.cali[,s] <= OBS.PREVALENCE.threshold,s] <- 0
          PA.OBS.PREVALENCE.eval[PA.OBS.PREVALENCE.eval[,s] >= OBS.PREVALENCE.threshold,s] <- 1
          PA.OBS.PREVALENCE.eval[PA.OBS.PREVALENCE.eval[,s] <= OBS.PREVALENCE.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/OBS.PREVALENCE", sep=""), recursive=TRUE)
      save(PA.OBS.PREVALENCE.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/OBS.PREVALENCE/PA.OBS.PREVALENCE.cali_",run,".RData", sep=""))
      save(PA.OBS.PREVALENCE.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/OBS.PREVALENCE/PA.OBS.PREVALENCE.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.OBS.PREVALENCE.cali,
                            PA.eval = PA.OBS.PREVALENCE.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/OBS.PREVALENCE/OBS.PREVALENCE_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by Avg Probability
    #-------------------------------------------------
    if("AVG.PROBABILITY" %in% thresholds){
      PA.AVG.PROBABILITY.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.AVG.PROBABILITY.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.AVG.PROBABILITY.cali)[2]){
        if(sum(!is.na(PA.AVG.PROBABILITY.cali[,s]))==0){
          PA.AVG.PROBABILITY.cali[,s] <- NA
          PA.AVG.PROBABILITY.eval[,s] <- NA
        }else{
          AVG.PROBABILITY.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.AVG.PROBABILITY.cali)[1]),
                                                                                  t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                                  t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                                          threshold=101, opt.methods=8)[,2]
          PA.AVG.PROBABILITY.cali[PA.AVG.PROBABILITY.cali[,s] >= AVG.PROBABILITY.threshold,s] <- 1
          PA.AVG.PROBABILITY.cali[PA.AVG.PROBABILITY.cali[,s] <= AVG.PROBABILITY.threshold,s] <- 0
          PA.AVG.PROBABILITY.eval[PA.AVG.PROBABILITY.eval[,s] >= AVG.PROBABILITY.threshold,s] <- 1
          PA.AVG.PROBABILITY.eval[PA.AVG.PROBABILITY.eval[,s] <= AVG.PROBABILITY.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/AVG.PROBABILITY", sep=""), recursive=TRUE)
      save(PA.AVG.PROBABILITY.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/AVG.PROBABILITY/PA.AVG.PROBABILITY.cali_",run,".RData", sep=""))
      save(PA.AVG.PROBABILITY.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/AVG.PROBABILITY/PA.AVG.PROBABILITY.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.AVG.PROBABILITY.cali,
                            PA.eval = PA.AVG.PROBABILITY.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/AVG.PROBABILITY/AVG.PROBABILITY_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by MCE
    #-------------------------------------------------
    if("MCE" %in% thresholds){
      PA.MCE.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.MCE.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      for(s in 1:dim(PA.MCE.cali)[2]){
        if(sum(!is.na(PA.MCE.cali[,s]))==0){
          PA.MCE.cali[,s] <- NA
          PA.MCE.eval[,s] <- NA
        }else{
          MCE.threshold <- optimal.thresholds(DATA=na.omit(data.frame(unlist(dimnames(PA.MCE.cali)[1]),
                                                                      t(ccv.modeling.data$speciesData.calibration[,,run])[,s],
                                                                      t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])[,s]/1000)), 
                                              threshold=101, opt.methods=10, req.sens=(100-MCE)/100)[,2]
          PA.MCE.cali[PA.MCE.cali[,s] >= MCE.threshold,s] <- 1
          PA.MCE.cali[PA.MCE.cali[,s] <= MCE.threshold,s] <- 0
          PA.MCE.eval[PA.MCE.eval[,s] >= MCE.threshold,s] <- 1
          PA.MCE.eval[PA.MCE.eval[,s] <= MCE.threshold,s] <- 0
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/MCE", sep=""), recursive=TRUE)
      save(PA.MCE.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MCE/PA.MCE.cali_",run,".RData", sep=""))
      save(PA.MCE.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MCE/PA.MCE.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.MCE.cali,
                            PA.eval = PA.MCE.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/MCE/MCE_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by PS_SDM
    #-------------------------------------------------
    if("PS_SDM" %in% thresholds){
      PA.PS_SDM.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.PS_SDM.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      #SR calculation
      SR.cali <- rowSums(PA.PS_SDM.cali, na.rm = TRUE)
      SR.eval <- rowSums(PA.PS_SDM.eval, na.rm = TRUE)
      for(p in 1:dim(PA.PS_SDM.cali)[1]){
        if(round(SR.cali[p])==0){
          PA.PS_SDM.cali[p,] <- 0
        }else{
          pS_SDM.threshold <- sort(PA.PS_SDM.cali[p,], decreasing = TRUE)[round(SR.cali[p])]
          PA.PS_SDM.cali[p,PA.PS_SDM.cali[p,]>=as.numeric(pS_SDM.threshold)] <- 1
          PA.PS_SDM.cali[p,PA.PS_SDM.cali[p,]<as.numeric(pS_SDM.threshold)] <- 0    
        }
      }
      for(p in 1:dim(PA.PS_SDM.eval)[1]){
        if(round(SR.eval[p])==0){
          PA.PS_SDM.eval[p,] <- 0
        }else{
          pS_SDM.threshold <- sort(PA.PS_SDM.eval[p,], decreasing = TRUE)[round(SR.eval[p])]
          PA.PS_SDM.eval[p,PA.PS_SDM.eval[p,]>=as.numeric(pS_SDM.threshold)] <- 1
          PA.PS_SDM.eval[p,PA.PS_SDM.eval[p,]<as.numeric(pS_SDM.threshold)] <- 0    
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/PS_SDM", sep=""), recursive = TRUE)
      save(PA.PS_SDM.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/PS_SDM/PA.PS_SDM.cali_",run,".RData", sep=""))
      save(PA.PS_SDM.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/PS_SDM/PA.PS_SDM.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.PS_SDM.cali,
                            PA.eval = PA.PS_SDM.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/PS_SDM/PS_SDM_", sep=""),
                            run=run)
    }
    
    #-------------------------------------------------
    #Thresholding by MEM
    #-------------------------------------------------
    if("MEM" %in% thresholds){
      PA.MEM.cali <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,run])/1000
      PA.MEM.eval <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,run])/1000
      #SR calculation
      SR.cali <- MEM[which(ccv.modeling.data$DataSplitTable[,run])]
      SR.eval <- MEM[which(!ccv.modeling.data$DataSplitTable[,run])]
      for(p in 1:length(SR.cali)){
        if(round(SR.cali[p])==0){
          PA.MEM.cali[p,] <- 0
        }else{
          MEM.threshold <- sort(PA.MEM.cali[p,], decreasing = TRUE)[round(SR.cali[p])]
          PA.MEM.cali[p,PA.MEM.cali[p,]>=as.numeric(MEM.threshold)] <- 1
          PA.MEM.cali[p,PA.MEM.cali[p,]<=as.numeric(MEM.threshold)] <- 0    
        }
      }
      for(p in 1:length(SR.eval)){
        if(round(SR.eval[p])==0){
          PA.MEM.eval[p,] <- 0
        }else{
          MEM.threshold <- sort(PA.MEM.eval[p,], decreasing = TRUE)[round(SR.eval[p])]
          PA.MEM.eval[p,PA.MEM.eval[p,]>=as.numeric(MEM.threshold)] <- 1
          PA.MEM.eval[p,PA.MEM.eval[p,]<=as.numeric(MEM.threshold)] <- 0    
        }
      }
      dir.create(paste(ccv.modeling.data$modeling.id, "/Thresholding/MEM", sep=""), recursive = TRUE)
      save(PA.MEM.cali, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MEM/PA.MEM.cali_",run,".RData", sep=""))
      save(PA.MEM.eval, file=paste(ccv.modeling.data$modeling.id, "/Thresholding/MEM/PA.MEM.eval_",run,".RData", sep=""))
      community.compairison(sp.data.cali=t(ccv.modeling.data$speciesData.calibration[,,run]), 
                            sp.data.eval= t(ccv.modeling.data$speciesData.evaluation[,,run]),
                            PA.cali= PA.MEM.cali,
                            PA.eval = PA.MEM.eval,
                            community.metrics = community.metrics,
                            save.dir=paste(ccv.modeling.data$modeling.id, "/Thresholding/MEM/MEM_", sep=""),
                            run=run)
    }
  }
  
  ##############################################################################################################
  ##############################################################################################################
  ##############################################################################################################
  
  #THE MAIN PROGRAMM
  if(parallel){
    #Run the stuff
    sfInit(parallel=TRUE, cpus=cpus)
    sfLibrary('PresenceAbsence', character.only=TRUE)
    sfExport('community.metrics.calculation')
    sfExport('community.compairison')
    sfLapply(1:dim(ccv.modeling.data$DataSplitTable)[2], 
             community.thresholding,
             ccv.modeling.data=ccv.modeling.data, 
             thresholds=thresholds, 
             community.metrics=community.metrics, 
             fix.threshold=fix.threshold, 
             MCE=MCE, 
             MEM=MEM)
    sfStop( nostop=FALSE )
  }else{
    lapply(1:dim(ccv.modeling.data$DataSplitTable)[2], 
           community.thresholding,
           ccv.modeling.data=ccv.modeling.data, 
           thresholds=thresholds, 
           community.metrics=community.metrics, 
           fix.threshold=fix.threshold, 
           MCE=MCE, 
           MEM=MEM)
  }
  
  #Gathering all the results
  ccv.metrics.allsites.cali <- array(data=NA,
                                     dim=c(dim(ccv.modeling.data$speciesData.full)[1],
                                           length(thresholds),
                                           length(community.metrics),
                                           dim(ccv.modeling.data$speciesData.calibration)[3]),
                                     dimnames=list(unlist(dimnames(ccv.modeling.data$speciesData.full)[1]),
                                                   thresholds,
                                                   community.metrics,
                                                   unlist(dimnames(ccv.modeling.data$speciesData.calibration)[3])))
  
  ccv.metrics.allsites.eval <- array(data=NA,
                                     dim=c(dim(ccv.modeling.data$speciesData.full)[1],
                                           length(thresholds),
                                           length(community.metrics),
                                           dim(ccv.modeling.data$speciesData.evaluation)[3]),
                                     dimnames=list(unlist(dimnames(ccv.modeling.data$speciesData.full)[1]),
                                                   thresholds,
                                                   community.metrics,
                                                   unlist(dimnames(ccv.modeling.data$speciesData.evaluation)[3])))
  
  ccv.PA.allSites <- array(data=NA,
                           dim=c(dim(ccv.modeling.data$speciesData.full)[2],
                                 dim(ccv.modeling.data$speciesData.full)[1],
                                 length(thresholds),
                                 dim(ccv.modeling.data$speciesData.evaluation)[3]),
                           dimnames=list(unlist(dimnames(ccv.modeling.data$speciesData.full)[2]),
                                         unlist(dimnames(ccv.modeling.data$speciesData.full)[1]),
                                         thresholds,
                                         unlist(dimnames(ccv.modeling.data$speciesData.evaluation)[3])))
  
  
  #Looping through all the files
  for(th in thresholds){
    for(r in 1:dim(ccv.modeling.data$DataSplitTable)[2]){

      
      #Loading the files
      load(paste(ccv.modeling.data$modeling.id,"/Thresholding/",th,"/",th,"_community.metrics.cali_",r,".RData", sep=""))
      load(paste(ccv.modeling.data$modeling.id,"/Thresholding/",th,"/",th,"_community.metrics.eval_",r,".RData", sep=""))
      load(paste(ccv.modeling.data$modeling.id,"/Thresholding/",th,"/",th,"_community.metrics.cali_",r,".RData", sep=""))
      load(paste(ccv.modeling.data$modeling.id,"/Thresholding/",th,"/PA.",th,".cali_",r,".RData", sep=""))
      load(paste(ccv.modeling.data$modeling.id,"/Thresholding/",th,"/PA.",th,".eval_",r,".RData", sep=""))
      
      #Parsing the results
      ccv.metrics.allsites.cali[which(ccv.modeling.data$DataSplitTable[,r]),th,,r] <- eval(parse(text="community.metrics.cali.selected"))[1:length(which(ccv.modeling.data$DataSplitTable[,r])),]
      
      ccv.metrics.allsites.eval[which(!ccv.modeling.data$DataSplitTable[,r]),th,,r] <- eval(parse(text="community.metrics.eval.selected"))[1:length(which(!ccv.modeling.data$DataSplitTable[,r])),]
      
      ccv.PA.allSites[,which(ccv.modeling.data$DataSplitTable[,r]),th,r] <- eval(parse(text=paste("PA.",th,".cali", sep="")))[1:length(which(ccv.modeling.data$DataSplitTable[,r]))]
      ccv.PA.allSites[,which(!ccv.modeling.data$DataSplitTable[,r]),th,r] <- eval(parse(text=paste("PA.",th,".eval", sep="")))[1:length(which(!ccv.modeling.data$DataSplitTable[,r]))]
    }
  }
  #Save and return
  ccv.evaluationMetrics.bin <- list(DataSplitTable = ccv.modeling.data$DataSplitTable,
                                    CommunityEvaluationMetrics.CalibrationSites = ccv.metrics.allsites.cali,
                                    CommunityEvaluationMetrics.EvaluationSites = ccv.metrics.allsites.eval,
                                    PA.allSites = ccv.PA.allSites)
  save(ccv.evaluationMetrics.bin, file=paste(ccv.modeling.data$modeling.id,".ccv.evaluationMetrics.bin.RData", sep=""))
  
  return(ccv.evaluationMetrics.bin)
}




#####################################################
# d) ecospat.CCV.communityEvaluation.prob
#####################################################

#DESCRIPTION
# evaluate the community predictions based directly on the probabilites (i.e., threshold indpendent) (ecospat.CCV.communityEvaluation.prob)
# This function generates a number of community evaluation metrics directly based on the probability/habitat suitability
# returned by the individual models.


#FUNCTION'S ARGUMENTS
#ccv.modeling.data     an output from ecospat.CCV.modeling function
#community.metric      probabilistic community metrics to calculate ("SR.deviation","community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard")
#parallel              binary variable if parallel computing is allowed
#cpus                  if parallel true the number of cpus to use


#DETAILS


#VALUES
#DataSplitTable                                           a matrix with TRUE/FALSE for each model run (TRUE=Calibration point, FALSE=Evaluation point)
#CommunityEvaluationMetrics.CalibrationSites              a 3-dimensional array containing the community evaluation metrics for the calibartion sites of each run (NA means that the site was used for evaluation)
#CommunityEvaluationMetrics.EvaluationSites               a 3-dimensional array containing the community evaluation metrics for the evaluation sites of each run (NA means that the site was used for calibaration)



#AUTHORS
# Daniel Scherrer <daniel.j.a.scherrer@gmail.com>


#REFERENCES
# Scherrer, D., D'Amen, M., Mateo, M.R.G., Fernandes, R.F. & Guisan , A. (2018) 
#   How to best threshold and validate stacked species assemblages? Community optimisation might hold the answer. 
#   Methods in Ecology and Evolution, 9(10): 2155-2166.
# Scherrer, D., Mod, H.K., Guisan, A. (2019)
#   How to evaluate community predictions without thresholding?
#   Methods in Ecology and Evolution, in press


#SEE ALSO
# ecospat.CCV.modelling

ecospat.CCV.communityEvaluation.prob <- function(ccv.modeling.data, 
                                                 community.metrics=c('SR.deviation','community.AUC','Max.Sorensen','Max.Jaccard','probabilistic.Sorensen','probabilistic.Jaccard'), 
                                                 parallel = FALSE, 
                                                 cpus = 4){
  
  
  #Loading the packages needed (they should all be installed by ecospat library)
  #require(PresenceAbsence)
  #require(poibin)
  #require(snowfall)
  
  #Checking all the input data
  stopifnot(names(ccv.modeling.data)==c("modeling.id",
                                        "output.files",
                                        "speciesData.calibration",
                                        "speciesData.evaluation",
                                        "speciesData.full",
                                        "DataSplitTable",
                                        "singleSpecies.ensembleEvaluationScore",
                                        "singleSpecies.ensembleVariableImportance",
                                        "singleSpecies.calibrationSites.ensemblePredictions",
                                        "singleSpecies.evaluationSites.ensemblePredictions",
                                        "allSites.averagePredictions.cali",
                                        "allSites.averagePredictions.eval"))
  
  stopifnot(community.metrics %in% c("SR.deviation","community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard"))
  
  #############################################################################################################
  #SUBFUNCTIONS################################################################################################
  #############################################################################################################
  
  #Community SR ##################################################################
  
  #For the SR.probability
  SR.prob <- function(data){
    Sj <- as.numeric(data[1])
    pjk <- as.numeric(data[-1][!is.na(data[-1])])
    return(dpoibin(kk=Sj, pp=pjk))
  }
  
  #For SR mean and sd
  SR.mean.sd <- function(data){
    data <- data[!is.na(data)]
    SR.mean <- sum(data[-1])
    SR.dev <- SR.mean - data[[1]]
    SR.sd <- sqrt(sum((1-data[-1])*data[-1]))
    if(SR.dev >= 0){
      SR.prob <- ppoibin(data[[1]], data[-1])
    }else{
      SR.prob <- 1-ppoibin(data[[1]]-1, data[-1])
    }
    return(unlist(c(SR.mean=SR.mean,SR.dev=SR.dev,SR.sd=SR.sd, SR.prob=SR.prob)))
  }
  
  #Community Composition ##################################################################
  
  #Community AUC
  Community.AUC <- function(data){
    obs.data <- as.numeric(data[1:(length(data)/2)])
    pred.data <- as.numeric(data[((length(data)/2)+1):length(data)])
    obs.data <- obs.data[!is.na(pred.data)]
    pred.data <- pred.data[!is.na(pred.data)]
    if(sum(is.na(obs.data))==length(obs.data) & sum(is.na(pred.data)==length(pred.data))){
      return(NA)
    }else{
      if(sum(obs.data)==0 | sum(obs.data)==length(obs.data)){
        return(1)
      }else{
        auc.return <- unlist(auc(DATA=data.frame(id=1:length(obs.data),
                                                 obs=obs.data,
                                                 pred=pred.data), na.rm = TRUE))[1]
        return(auc.return)
      }
    }
  }
  
  #For the community.probability
  composition.prob <- function(data){
    obs.data <- data[1:(length(data)/2)]
    pred.data <- data[((length(data)/2)+1):length(data)]
    obs.data <- obs.data[!is.na(pred.data)]
    pred.data <- pred.data[!is.na(pred.data)]
    if(sum(obs.data==1)>0 & sum(obs.data==0)>0){
      prob.list <- c(pred.data[which(obs.data==1)],1-pred.data[which(obs.data==0)])
    }
    if(sum(obs.data==1)>0 & sum(obs.data==0)==0){
      prob.list <- pred.data[which(obs.data==1)]
    }
    if(sum(obs.data==1)==0 & sum(obs.data==0)>0){
      prob.list <- 1-pred.data[which(obs.data==0)]
    }
    return(prod(prob.list))
  }
  
  #For the Max.Sorensen
  MaxSorensen <- function(data){
    obs.data <- as.numeric(data[1:(length(data)/2)])
    pred.data <- as.numeric(data[((length(data)/2)+1):length(data)])
    temp.Sorensen <- rep(NA,101)
    th <- seq(0,1,0.01)
    for(i in 1:101){
      pred.temp <- pred.data
      pred.temp[pred.temp>=th[i]] <- 1
      pred.temp[pred.temp<th[i]] <- 0
      errors <- 2*pred.temp+obs.data
      a <- length(which(errors == 3)) #True Positives
      b <- length(which(errors == 2)) #False Positives
      c <- length(which(errors == 1)) #False Negatives
      if(a==0 & b==0 & c==0){
        Sorensen <- 1
      }else{
        Sorensen <- round((2 * a)/(2 * a + b + c), digits=3)
      }
      temp.Sorensen[i] <- Sorensen
    }
    return(max(temp.Sorensen))
  }
  
  #For the Max.Jaccard
  MaxJaccard <- function(data){
    obs.data <- as.numeric(data[1:(length(data)/2)])
    pred.data <- as.numeric(data[((length(data)/2)+1):length(data)])
    temp.Jaccard <- rep(NA,101)
    th <- seq(0,1,0.01)
    for(i in 1:101){
      pred.temp <- pred.data
      pred.temp[pred.temp>=th[i]] <- 1
      pred.temp[pred.temp<th[i]] <- 0
      errors <- 2*pred.temp+obs.data
      a <- length(which(errors == 3)) #True Positives
      b <- length(which(errors == 2)) #False Positives
      c <- length(which(errors == 1)) #False Negatives
      if(a==0 & b==0 & c==0){
       Jaccard <- 1
      }else{
        Jaccard <- round((a)/(a + b + c), digits=3)
      }
      temp.Jaccard[i] <- Jaccard
    }
    return(max(temp.Jaccard))
  }
  
  #For probabilistic Sorensen
  probabilisticSorensen <- function(data){
    temp.df <- data.frame(obs=as.numeric(data[1:(length(data)/2)]),pred=as.numeric(data[((length(data)/2)+1):length(data)]))
    temp.df <- temp.df[order(-temp.df$pred),]
    AnB <- 2* sum(temp.df$pred[temp.df$obs==1])
    AuB <- sum(temp.df$pred[temp.df$pred>=min(temp.df$pred[temp.df$obs==1])]) + sum(temp.df$pred[temp.df$obs==1])
    return(AnB/AuB)
  }
  
  #For probabilistic Jaccard
  probabilisticJaccard <- function(data){
    temp.df <- data.frame(obs=as.numeric(data[1:(length(data)/2)]),pred=as.numeric(data[((length(data)/2)+1):length(data)]))
    temp.df <- temp.df[order(-temp.df$pred),]
    AnB <- sum(temp.df$pred[temp.df$obs==1])
    AuB <- sum(temp.df$pred[temp.df$pred>=min(temp.df$pred[temp.df$obs==1])])
    return(AnB/AuB)
  }
  
  #################################################################################################
  #The main sub-function to calculate the metrics##################################################
  #################################################################################################
  prob.community.metics <- function(obs, pred, metrics){
    
    #Some basic arrangements
    obs <- obs[,order(colnames(obs))]
    pred <- pred[,order(colnames(pred))]
    
    #Test the input data
    try(if(!identical(dim(obs),dim(pred))){stop("Dimensions of obs and pred differ")})
    try(if(!identical(colnames(obs), colnames(pred))){stop("Columnames of obs and pred differ, make sure the species are matching")})
    
    #Making the Null models###############################################################
    #Null model 1. Only the species pool and number of sites is known
    Null.pred.05 <- pred
    Null.pred.05[,] <- 0.5
    #Null model 2. Species pool, number of sites and average SR is known
    Null.pred.average.SR <- pred
    Null.pred.average.SR[,] <- mean(rowSums(obs))/dim(obs)[2]
    #Reflecting the observed prevalence
    Null.pred.prevalence <- pred
    Null.pred.prevalence[,] <- rep(colSums(obs)/dim(obs)[1], each=dim(obs)[1])
    
    #Calculating the SR metrics
    if("SR.deviation" %in% metrics){
      #The observed SR
      SR.obs <- rowSums(obs)
      #The three NULL predictions
      SR.Null.pred.05 <- apply(data.frame(SR.obs=SR.obs, Null.pred.05),1, SR.prob)
      SR.Null.pred.average.SR <- apply(data.frame(SR.obs=SR.obs, Null.pred.average.SR),1, SR.prob)
      SR.Null.pred.prevalence <- apply(data.frame(SR.obs=SR.obs, Null.pred.prevalence),1, SR.prob)
      #The actual probability
      SR.pred <- apply(data.frame(SR.obs=SR.obs, pred),1,SR.prob)
      #The mean and sd of the SR
      SR.stat <- data.frame(t(apply(data.frame(SR.obs=SR.obs, pred),1,SR.mean.sd)))
      #The improvement over the two NULL models
      SR.results <- signif(data.frame(SR.obs, SR.stat, SR.imp.05=SR.pred/SR.Null.pred.05, SR.imp.average.SR=SR.pred/SR.Null.pred.average.SR, SR.imp.prevalence=SR.pred/SR.Null.pred.prevalence),3)
    }
    
    #Calculating the composition metrics
    if(length(intersect(metrics,c("community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard")))>0){
      
      #The two NULL predictions
      composition.Null.pred.05 <- rep(0.5^dim(Null.pred.05)[2], dim(Null.pred.05)[1])
      composition.Null.pred.average.SR <- apply(data.frame(obs, Null.pred.average.SR),1, composition.prob)
      composition.Null.pred.prevalence <- apply(data.frame(obs, Null.pred.prevalence),1, composition.prob)
      
      #The actual probability
      composition.pred <- apply(data.frame(obs, pred),1, composition.prob)
      composition.results <- signif(data.frame(composition.imp.05 = composition.pred/composition.Null.pred.05, composition.imp.average.SR = composition.pred/composition.Null.pred.average.SR, composition.imp.prevalence = composition.pred/composition.Null.pred.prevalence),3)
      
      #The probabilistic.Sorensen
      if("probabilistic.Sorensen" %in% metrics){
        Sorensen.stat <- data.frame(t(apply(data.frame(obs, pred),1, probabilisticSorensen)))
        composition.results <- signif(data.frame(probabilistic.Sorensen=unlist(Sorensen.stat), composition.results),3)
      }
      
      #The probabilistic.Jaccard
      if("probabilistic.Jaccard" %in% metrics){
        Jaccard.stat <- data.frame(t(apply(data.frame(obs, pred),1, probabilisticJaccard)))
        composition.results <- signif(data.frame(probabilistic.Jaccard=unlist(Jaccard.stat), composition.results),3)
      }
      
      #The Max.Sorensen
      if("Max.Sorensen" %in% metrics){
        Sorensen.stat <- data.frame(t(apply(data.frame(obs, pred),1, MaxSorensen)))
        composition.results <- signif(data.frame(Max.Sorensen=unlist(Sorensen.stat), composition.results),3)
      }
      
      #The Max.Jaccard
      if("Max.Jaccard" %in% metrics){
        Jaccard.stat <- data.frame(t(apply(data.frame(obs, pred),1, MaxJaccard)))
        composition.results <- signif(data.frame(Max.Jaccard=unlist(Jaccard.stat), composition.results),3)
      }
      
      #The mean and sd AUC
      if("community.AUC" %in% metrics){
        AUC.stat <- apply(data.frame(obs, pred),1, Community.AUC)
        composition.results <- signif(data.frame(Community.AUC=AUC.stat, composition.results),3)
      }
    }
    
    #Returning the results
    if("SR.deviation" %in% metrics & length(intersect(metrics,c("community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard")))==0){
      return(SR.results)
    }
    if(!("SR.deviation" %in% metrics) & length(intersect(metrics,c("community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard")))>0){
      return(composition.results)
    }
    if("SR.deviation" %in% metrics & length(intersect(metrics,c("community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard")))>0){
      return(data.frame(SR.results,composition.results))
    }
  }
  
  #################################################################################################
  #THE MAIN PROGRAMM###############################################################################
  #################################################################################################
  #Calculation of array width (nb of measurements)
  nb.mes <- 0
  if("SR.deviation" %in% community.metrics){
    nb.mes <- nb.mes+8
  }
  if(length(intersect(community.metrics,c("community.AUC","Max.Sorensen","Max.Jaccard","probabilistic.Sorensen","probabilistic.Jaccard")))>0){
    nb.mes <- nb.mes+3  
  }
  if("community.AUC" %in% community.metrics){
    nb.mes <- nb.mes+1 
  }
  if("Max.Sorensen" %in% community.metrics){
    nb.mes <- nb.mes+1 
  }
  if("Max.Jaccard" %in% community.metrics){
    nb.mes <- nb.mes+1 
  }
  if("probabilistic.Sorensen" %in% community.metrics){
    nb.mes <- nb.mes+1 
  }
  if("probabilistic.Jaccard" %in% community.metrics){
    nb.mes <- nb.mes+1 
  }

  
  #Array to save data
  ccv.cali <- array(data=NA, dim=c(dim(ccv.modeling.data$speciesData.calibration)[2],nb.mes, dim(ccv.modeling.data$speciesData.calibration)[3]))
  ccv.eval <- array(data=NA, dim=c(dim(ccv.modeling.data$speciesData.evaluation)[2],nb.mes, dim(ccv.modeling.data$speciesData.evaluation)[3]))
  
  #Making the computations
  if(parallel){
    sfInit(parallel=TRUE, cpus=cpus)
    sfExport("SR.mean.sd", "SR.prob","prob.community.metics", "composition.prob","Community.AUC", "MaxJaccard", "MaxSorensen", "probabilisticJaccard", "probabilisticSorensen")
    sfExport("ccv.modeling.data", "community.metrics")
    sfLibrary("poibin", character.only=TRUE )
    sfLibrary("PresenceAbsence", character.only=TRUE )
    
    #Calibration data
    temp <- sfLapply(1:dim(ccv.modeling.data$speciesData.calibration)[3], 
                     function(x){
                       obs.temp <- t(ccv.modeling.data$speciesData.calibration[,,x])[rowSums(is.na(t(ccv.modeling.data$speciesData.calibration[,,x])))!=ncol(t(ccv.modeling.data$speciesData.calibration[,,x])),]
                       pred.temp <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,x])[rowSums(is.na(t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,x])))!=ncol(t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,x])),]/1000
                       stopifnot(dim(obs.temp)==dim(pred.temp))
                       prob.community.metics(obs=obs.temp, 
                                             pred=pred.temp, 
                                             metrics=community.metrics)
                     })
    for(i in 1:dim(ccv.modeling.data$speciesData.calibration)[3]){
      ccv.cali[1:dim(temp[[i]])[1],,i] <- unlist(temp[[i]])
    }
    dimnames(ccv.cali) <- list(dimnames(ccv.modeling.data$speciesData.calibration)[[2]], unlist(dimnames(temp[[1]])[2]), dimnames(ccv.modeling.data$speciesData.calibration)[[3]])
    
    #Evaluation data
    temp <- sfLapply(1:dim(ccv.modeling.data$speciesData.evaluation)[3], 
                     function(x){
                       obs.temp <- t(ccv.modeling.data$speciesData.evaluation[,,x])[rowSums(is.na(t(ccv.modeling.data$speciesData.evaluation[,,x])))!=ncol(t(ccv.modeling.data$speciesData.evaluation[,,x])),]
                       pred.temp <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,x])[rowSums(is.na(t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,x])))!=ncol(t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,x])),]/1000
                       stopifnot(dim(obs.temp)==dim(pred.temp))
                       prob.community.metics(obs=obs.temp, 
                                             pred=pred.temp, 
                                             metrics=community.metrics)
                     })
    for(i in 1:dim(ccv.modeling.data$speciesData.evaluation)[3]){
      ccv.eval[1:dim(temp[[i]])[1],,i] <- unlist(temp[[i]])
    }
    dimnames(ccv.eval) <- list(dimnames(ccv.modeling.data$speciesData.evaluation)[[2]], unlist(dimnames(temp[[1]])[2]), dimnames(ccv.modeling.data$speciesData.evaluation)[[3]])
    sfStop( nostop=FALSE )
  }else{
    
    #Calibration data
    temp <- lapply(1:dim(ccv.modeling.data$speciesData.calibration)[3], 
                   function(x){
                     obs.temp <- t(ccv.modeling.data$speciesData.calibration[,,x])[rowSums(is.na(t(ccv.modeling.data$speciesData.calibration[,,x])))!=ncol(t(ccv.modeling.data$speciesData.calibration[,,x])),]
                     pred.temp <- t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,x])[rowSums(is.na(t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,x])))!=ncol(t(ccv.modeling.data$singleSpecies.calibrationSites.ensemblePredictions[,,x])),]/1000
                     stopifnot(dim(obs.temp)==dim(pred.temp))
                     prob.community.metics(obs=obs.temp, 
                                           pred=pred.temp, 
                                           metrics=community.metrics)
                   })
    for(i in 1:dim(ccv.modeling.data$speciesData.calibration)[3]){
      ccv.cali[1:dim(temp[[i]])[1],,i] <- unlist(temp[[i]])
    }
    dimnames(ccv.cali) <- list(dimnames(ccv.modeling.data$speciesData.calibration)[[2]], unlist(dimnames(temp[[1]])[2]), dimnames(ccv.modeling.data$speciesData.calibration)[[3]])
    
    #Evaluation data
    temp <- lapply(1:dim(ccv.modeling.data$speciesData.evaluation)[3], 
                   function(x){
                     obs.temp <- t(ccv.modeling.data$speciesData.evaluation[,,x])[rowSums(is.na(t(ccv.modeling.data$speciesData.evaluation[,,x])))!=ncol(t(ccv.modeling.data$speciesData.evaluation[,,x])),]
                     pred.temp <- t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,x])[rowSums(is.na(t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,x])))!=ncol(t(ccv.modeling.data$singleSpecies.evaluationSites.ensemblePredictions[,,x])),]/1000
                     stopifnot(dim(obs.temp)==dim(pred.temp))
                     prob.community.metics(obs=obs.temp, 
                                           pred=pred.temp, 
                                           metrics=community.metrics)
                   })
    for(i in 1:dim(ccv.modeling.data$speciesData.evaluation)[3]){
      ccv.eval[1:dim(temp[[i]])[1],,i] <- unlist(temp[[i]])
    }
    dimnames(ccv.eval) <- list(dimnames(ccv.modeling.data$speciesData.evaluation)[[2]], unlist(dimnames(temp[[1]])[2]), dimnames(ccv.modeling.data$speciesData.evaluation)[[3]])
  }
  
  #For all sites
  CommunityEvaluationMetrics.CalibrationSites <- array(data=NA, 
                             dim=c(dim(ccv.modeling.data$speciesData.full)[1],nb.mes, dim(ccv.modeling.data$speciesData.calibration)[3]),
                             dimnames=list(unlist(dimnames(ccv.modeling.data$speciesData.full)[1]), unlist(dimnames(ccv.cali)[2]), unlist(dimnames(ccv.cali)[3])))
  CommunityEvaluationMetrics.EvaluationSites <- array(data=NA, 
                             dim=c(dim(ccv.modeling.data$speciesData.full)[1],nb.mes, dim(ccv.modeling.data$speciesData.calibration)[3]),
                             dimnames=list(unlist(dimnames(ccv.modeling.data$speciesData.full)[1]), unlist(dimnames(ccv.cali)[2]), unlist(dimnames(ccv.cali)[3])))
  
  for(r in 1:dim(ccv.modeling.data$speciesData.calibration)[3]){
    CommunityEvaluationMetrics.CalibrationSites[which(ccv.modeling.data$DataSplitTable[,r]),,r] <- ccv.cali[1:sum(ccv.modeling.data$DataSplitTable[,r]),,r]
    CommunityEvaluationMetrics.EvaluationSites[which(!ccv.modeling.data$DataSplitTable[,r]),,r] <- ccv.eval[1:sum(!ccv.modeling.data$DataSplitTable[,r]),,r]
  }
  
  #Save and return
  ccv.evaluationMetrics.prob <- list(DataSplitTable = ccv.modeling.data$DataSplitTable,
                                     CommunityEvaluationMetrics.CalibrationSites=CommunityEvaluationMetrics.CalibrationSites,
                                     CommunityEvaluationMetrics.EvaluationSites=CommunityEvaluationMetrics.EvaluationSites)
  save(ccv.evaluationMetrics.prob, file=paste(ccv.modeling.data$modeling.id,".ccv.evaluationMetrics.prob.RData", sep=""))
  
  return(ccv.evaluationMetrics.prob)
}
