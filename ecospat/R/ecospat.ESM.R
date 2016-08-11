#####################################################################################
# Commented functions for calculating  "Ensemble of Small Models" (ESM) in R using
# the Biomod2 package
#
# Version 11 [20/05/2015]
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
## models:              vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF','MAXENT.Phillips', "MAXENT.Tsuruoka" (same as in Biomod)
## modeling.id:	        character, the ID (=name) of modeling procedure. A random number by default.
## models.options:      BIOMOD.models.options object returned by BIOMOD_ModelingOptions (same as in Biomod)
## which.biva:          integer. which bivariate combinations should be used for modeling? Default: all 
## weighting.score:     evaluation score used to weight single models to build ensembles: 'AUC', 'SomersD' (2xAUC-1), 'Kappa', 'TSS' or 'Boyce'
## ESM_Projection:      logical. set to FALSE, if no projections should be calculated. Only Evaluation scores based on test and train data will be returned. (Default: TRUE)
## new.env:             A set of explanatory variables onto which models will be projected . It could be a data.frame, a matrix, or a rasterStack object. Make sure the column names (data.frame or matrix) or layer Names (rasterStack) perfectly match with the names of variables used to build the models in the previous steps.
## parallel:            logical. If TRUE, the parallel computing is enabled (highly recommended)
## cleanup:             numeric. Calls removeTmpFiles() to delete all files from rasterOptions()$tmpdir which are older than the given time (in hours). This is might be necessary to prevent running over quota. No cleanup is used by default.  


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
# Mirko Di Febbraro <mirkodifebbraro@gmail.com>
# Frank Breiner
# Olivier Broennimann


## References
#Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modeling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
# Breiner F.T., Guisan A., Bergamini A., Nobis M.P. (2015):
#  Overcoming limitations of modeling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218. doi: 10.1111/2041-210X.12403

##See Also
#ecospat.ESM.EnsembleModeling; ecospat.ESM.MergeModels

ecospat.ESM.Modeling <- function(data, NbRunEval=NULL, DataSplit, DataSplitTable=NULL, weighting.score,
                    models, modeling.id=as.character(format(Sys.time(), '%s')), models.options, 
                    which.biva=NULL, parallel, cleanup=F)
{

  #if(!require(biomod2)){stop("biomod2 package required!")}
  #if(!require(gtools)){stop("gtools package required!")}
  
  ## detach package GAM
  #if("GAM"%in%models){ detach(package:ecospat,force=T);detach(package:gam,force=T);unloadNamespace("ecospat")}
  
  if(!weighting.score %in% c("AUC","TSS","Boyce","Kappa","SomersD")){stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")}
  if(data@has.data.eval){stop("Evaluation with independant data is not supported yet!")}
  if("PA" %in% slotNames(data)){if(ncol(data@PA)>1){stop("It is not possible to use more than one Pseudo Absences dataset")}}
  
  if ("MAXENT.Phillips" %in% models) {
    if (!file.exists(paste(myBiomodOption@MAXENT.Phillips$path_to_maxent.jar,"maxent.jar",sep="/"))) {
      stop("maxent.jar file not found!")}}
  
  if(weighting.score =='AUC' | weighting.score =='SomersD'){
  models.eval.meth <- 'ROC'}    
  if(weighting.score =='Kappa'){
  models.eval.meth <- 'KAPPA'}
  if(weighting.score =='TSS'){
  models.eval.meth <- 'TSS'}
  if(weighting.score =='Boyce'){
  models.eval.meth <- 'ROC'}    
  
  models <- sort(models)
  
  iniwd<-getwd()
  dir.create(paste("./ESM.BIOMOD.output",data@sp.name,sep="_"))
  newwd <- paste(getwd(),"/ESM.BIOMOD.output_",data@sp.name,sep="")
  setwd(newwd)
 
  combinations<-combn(colnames(data@data.env.var), 2)
  if(is.null(which.biva)){ which.biva <- 1:ncol(combinations)}

  mydata<-data

  # produce calib.lines and keep DataSplitTable constant for each BiVa Model
  if(is.null(DataSplitTable)){
  mod.prep.dat <- .Models.prepare.data(mydata, NbRunEval, DataSplit, Yweights=NULL, Prevalence=0.5, do.full.models=T)
  if(length(dim(mod.prep.dat[[1]]$calibLines))==3){
  	calib.lines <- mod.prep.dat[[1]]$calibLines[,,1]
  }
  if(length(dim(mod.prep.dat[[1]]$calibLines))==2){
  	calib.lines <- mod.prep.dat[[1]]$calibLines
  }
  rm(mod.prep.dat)
  }else{calib.lines <- DataSplitTable}
  
	if(is.null(NbRunEval)){
		if(ncol(calib.lines>1)){
			if(sum(!calib.lines[,ncol(calib.lines)])==0){
				NbRunEval <-  ncol(calib.lines)-1
			}else{
				NbRunEval <-  ncol(calib.lines)
			}	
		}else{
			if(sum(!calib.lines)==0){
				NbRunEval <-  ncol(calib.lines)-1
			}else{
				NbRunEval <-  ncol(calib.lines)
			}
		}
	}


  mymodels<-list()
           
  if(parallel==F)
  {
    for(k in which.biva)
    {

      mydata@data.env.var <- data@data.env.var[,colnames(data@data.env.var)%in%combinations[,k]]
      mydata@sp.name <- paste("ESM.BIOMOD", k, sep=".")
      
      mymodels[[k]]<-BIOMOD_Modeling(data=mydata, 
                               models=models,
                               models.options=models.options, 
                               models.eval.meth=models.eval.meth,
                               DataSplitTable=calib.lines,
                               Prevalence=0.5,
                               rescal.all.models=T,
                               do.full.models=T,
                               VarImport=0,
                               modeling.id=modeling.id)
                                     
      if(cleanup!=F){
        removeTmpFiles(h=cleanup)
      }
    
    }
  }
  
  if(parallel==T)
  {
        
    mymodels<-foreach(k=which.biva, .packages="biomod2") %dopar%
    {
      setwd(newwd)
      mydata@data.env.var <- data@data.env.var[,colnames(data@data.env.var)%in%combinations[,k]]
      mydata@sp.name <- paste("ESM.BIOMOD", k, sep=".")
     
      if(cleanup!=F){
        removeTmpFiles(h=cleanup)
      }
      
      BIOMOD_Modeling(data=mydata, 
                               models=models,
                               models.options=models.options, 
                               models.eval.meth=models.eval.meth,
                               DataSplitTable=calib.lines,
                               Prevalence=0.5,
                               rescal.all.models=T,
                               do.full.models=T,
                               VarImport=0,
				                       modeling.id=modeling.id)
  
    }
  }

  output <- list(modeling.id = modeling.id,
        models. = grep(modeling.id,mixedsort(list.files(getwd(), "models.out", recursive=T, full.names=T)),value=T),
        models = models,
        calib.lines = calib.lines,
        NbRunEval = NbRunEval,
        data = data,
        wd = getwd(),
        which.biva = which.biva,
        mymodels=mymodels)
        
    save(output, file=paste("ESM_Modeling..models",modeling.id,"out",sep="."))
## attach package gam
#if("GAM"%in%models){ library(gam);library(ecospat)} 
  
    setwd(iniwd)
  return(output)
}


############################
# b)  ecospat.ESM.Projection    projects simple bivariate models on new.env.
############################

## FUNCTION'S ARGUMENTS
## ESM.modeling.output:   BIOMOD.formated.data object returned by BIOMOD_FormatingData
## new.env:             A set of explanatory variables onto which models will be projected . It could be a data.frame, a matrix, or a rasterStack object. Make sure the column names (data.frame or matrix) or layer Names (rasterStack) perfectly match with the names of variables used to build the models in the previous steps.
## parallel:            logical. If TRUE, the parallel computing is enabled

## Note: The name of new.env must be a regular expression (see ?regex)

ecospat.ESM.Projection <- function(ESM.modeling.output, new.env,parallel=F){
  
  iniwd<-getwd()
  setwd(ESM.modeling.output$wd)
  models<-ESM.modeling.output$models
  models.<-ESM.modeling.output$models.
  mymodels<-ESM.modeling.output$mymodels
  combinations<-combn(colnames(ESM.modeling.output$data@data.env.var), 2)
  which.biva <- ESM.modeling.output$which.biva
  NbRunEval<-ESM.modeling.output$NbRunEval
  modeling.id<-ESM.modeling.output$modeling.id
  name.env<-deparse(substitute(new.env))
## detach package GAM
#if("GAM"%in%models){ detach(package:ecospat,force=T);detach(package:gam,force=T);unloadNamespace("ecospat")}


  if(parallel==F)
  {
  for (k in 1:length(mymodels)){
    mymodel<-mymodels[[k]]
    if(is.data.frame(new.env)){
      BIOMOD_Projection(modeling.output=mymodel,
                                new.env=new.env[,colnames(new.env)%in%combinations[,k]],
                                proj.name=paste(name.env,"ESM.BIOMOD", k,mymodel@modeling.id, sep="."),
                                selected.models = c(grep("Full",mymodel@models.computed,value=T),
                                                    grep(paste("RUN",NbRunEval+1,sep=""),mymodel@models.computed,value=T)), # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="") 
                                do.stack=F,
                                build.clamping.mask = F)
    }
    
    if(class(new.env)=="RasterStack"){
      BIOMOD_Projection(modeling.output=mymodel,
                                new.env=new.env[[which(names(new.env)%in%combinations[,k])]],
                                proj.name=paste(name.env,"ESM.BIOMOD", k,mymodel@modeling.id, sep="."),
                                selected.models = c(grep("Full",mymodel@models.computed,value=T),
                                                    grep(paste("RUN",NbRunEval+1,sep=""),mymodel@models.computed,value=T)), # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="") 
                                do.stack=T,
                                build.clamping.mask = F)
    }
  }
  }
  if(parallel==T){
      foreach(k=1:length(mymodels), .packages="biomod2") %dopar%{
      mymodel<-mymodels[[k]]
      if(is.data.frame(new.env)){
        BIOMOD_Projection(modeling.output=mymodel,
                                  new.env=new.env[,colnames(new.env)%in%combinations[,k]],
                                  proj.name=paste(name.env,"ESM.BIOMOD", k,mymodel@modeling.id, sep="."),
                                  selected.models = c(grep("Full",mymodel@models.computed,value=T),
                                                      grep(paste("RUN",NbRunEval+1,sep=""),mymodel@models.computed,value=T)), # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="") 
                                  do.stack=F,
                                  build.clamping.mask = F)
      }
      
      if(class(new.env)=="RasterStack"){
        BIOMOD_Projection(modeling.output=mymodel,
                                  new.env=new.env[[which(names(new.env)%in%combinations[,k])]],
                                  proj.name=paste(name.env,"ESM.BIOMOD", k,mymodel@modeling.id, sep="."),
                                  selected.models = c(grep("Full",mymodel@models.computed,value=T),
                                                      grep(paste("RUN",NbRunEval+1,sep=""),mymodel@models.computed,value=T)), # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="") 
                                  do.stack=T,
                                  build.clamping.mask = F)
      }
    }
    
  }
  
    
  
  removeTmpFiles(h=T)
  
  pred.biva<-grep(modeling.id,mixedsort(list.files(getwd(), paste("proj_",name.env,sep=""), recursive=T, full.names=T)) ,value=T)        
  
  output <- list(modeling.id = modeling.id,
                 models. = grep(modeling.id,mixedsort(list.files(getwd(), "models.out", recursive=T, full.names=T)),value=T),
                 models = models,
                 pred.biva = grep(modeling.id,mixedsort(list.files(getwd(), paste("proj_",name.env,sep=""), recursive=T, full.names=T)) ,value=T),
                 NbRunEval = NbRunEval,
                 name.env = name.env,
                 new.env.raster = class(new.env)=="RasterStack", 
                 wd = getwd(),
                 which.biva = which.biva)
  
  save(output, file=paste("ESM_Projections",modeling.id,"out",sep="."))
## attach package gam
#if("GAM"%in%models){ library(gam);library(ecospat)}

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
# ESM.predictions:  Returns the projections of ESMs for the selected single models and their ensemble

## Authors:
# Frank Breiner

## References
#Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
#Breiner, F. B., A. Guisan, A. Bergamini, M. P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, submitted.

## See also
# ecospat.ESM.Modeling; ecospat.ESM.MergeModels 

ecospat.ESM.EnsembleModeling <- function(ESM.modeling.output, weighting.score,threshold=NULL,
                    models)
{
 
  #if(!require(biomod2)){stop("biomod2 package required!")}
  #if(!require(PresenceAbsence)){stop("PresenceAbsence package required!")}
  if(!weighting.score %in% c("AUC","TSS","Boyce","Kappa","SomersD")){stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")}
  #if(!require(gtools)){stop("gtools package required!")}  
  #if(!require(ecospat)){stop("ecospat package required!")} 
   
  iniwd<-getwd()
  setwd(ESM.modeling.output$wd)
  data <- ESM.modeling.output$data
  models. <- ESM.modeling.output$models.
  NbRunEval <-  ESM.modeling.output$NbRunEval
  models <-  ESM.modeling.output$models
  calib.lines <-  ESM.modeling.output$calib.lines
  ESM_Projection <- ESM.modeling.output$ESM_Projection
  new.env <-  ESM.modeling.output$new.env
  
  for(i in 1:length(models.))load(models.[i])
 
models.<-NULL
  for(n in 1:length(ESM.modeling.output$modeling.id)){
    models.<-c(models.,mixedsort(grep(ESM.modeling.output$modeling.id[n], grep(paste(".",ESM.modeling.output$modeling.id[n],".models.out",sep=""), ls(), value=T,fixed=T), value=T)))
    }
    
  mymodel<-list()

  for(i in 1:length(models.))mymodel[[i]]<-get(models.[i])

### build the weighting vector to average the bivariate predictions   
     weights<-unlist(lapply(mymodel, function(x){
        y<-x
        if(weighting.score =='Boyce'){
          z <-get_predictions(y,as.data.frame=T) 
          z <- z[,grep(paste("RUN",NbRunEval+1,sep=""),colnames(z),invert=T)]
          x<-get_evaluations(y)[,"Testing.data",,,]
          if(length(models)>1){    
            x[,]<-NA
            x<-x[,colnames(x)!="Full" & colnames(x)!=paste("RUN",NbRunEval+1,sep="")] # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="")          
          }else{ 
            x[]<-NA
            x<-x[names(x)!="Full" & names(x)!=paste("RUN",NbRunEval+1,sep="")] # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="")          
          }  
          for(n in 1:(length(models))){
            model <- models[n]  
            for(i in 1:NbRunEval){
              if(sum(is.na(z[,grep(paste("RUN",i,"_",model,sep=""),colnames(z))]))!=nrow(z)){
              if(length(models)>1){
                # if(packageVersion("ecospat")==1.0){
                #   x[rownames(x)==model,colnames(x)==paste("RUN",i,sep="")] <- ecospat.boyce(z[!calib.lines[,paste("_RUN",i,sep="")],grep(paste("RUN",i,"_",model,sep=""),colnames(z))], z[!calib.lines[,paste("_RUN",i,sep="")] & data@data.species==1, grep(paste("RUN",i,"_",model,sep=""),colnames(z))], PEplot = F)$Pearson.cor
                # }else{
                  x[rownames(x)==model,colnames(x)==paste("RUN",i,sep="")] <- ecospat.boyce(z[!calib.lines[,paste("_RUN",i,sep="")],grep(paste("RUN",i,"_",model,sep=""),colnames(z))], z[!calib.lines[,paste("_RUN",i,sep="")] & data@data.species==1, grep(paste("RUN",i,"_",model,sep=""),colnames(z))], PEplot = F)$Spearman.cor
                #}           
              }else{ 
                # if(packageVersion("ecospat")==1.0){
                #  x[names(x)==paste("RUN",i,sep="")] <- ecospat.boyce(z[!calib.lines[,paste("_RUN",i,sep="")],grep(paste("RUN",i,"_",model,sep=""),colnames(z))], z[!calib.lines[,paste("_RUN",i,sep="")] & data@data.species==1, grep(paste("RUN",i,"_",model,sep=""),colnames(z))], PEplot = F)$Pearson.cor
                # }else{
                  x[names(x)==paste("RUN",i,sep="")] <- ecospat.boyce(z[!calib.lines[,paste("_RUN",i,sep="")],grep(paste("RUN",i,"_",model,sep=""),colnames(z))], z[!calib.lines[,paste("_RUN",i,sep="")] & data@data.species==1, grep(paste("RUN",i,"_",model,sep=""),colnames(z))], PEplot = F)$Spearman.cor
                #}
              }    
            }}
          }
          if(length(models)>1){          
            x<-round(apply(x,1,mean,na.rm=T),4)
          }else{ 
            x<-round(mean(x,na.rm=T),4)  
          }
        }else{
        x<-get_evaluations(y)[,"Testing.data",,,]
        	if(length(models)>1){
        	  for(row in models){
        	   if(ncol(calib.lines)==NbRunEval){ ## otherwise error in the next if bracket when no full model was calibrated
        	    if(is.na(x[row,NbRunEval+1])){   ## if NbRunEval+1 is na, the model using full data failed
        	      x[row,]<-NA
        	    }}
        	  }             
          		x<-x[,colnames(x)!="Full" & colnames(x)!=paste("RUN",NbRunEval+1,sep="")] # if DataSplitTable is provided to BIOMOD_Modeling, Full models are named: paste("RUN",NbRunEval+1,sep="") 
          	if(NbRunEval>1){
              x<-round(apply(x,1,mean,na.rm=T),4)
          	}
        	}else{
        	 if(ncol(calib.lines)==NbRunEval){ ## otherwise error in the next if bracket when no full model was calibrated
        	  if(is.na(x[NbRunEval+1])){       ## if NbRunEval+1 is na, the model using full data failed
        	    x<-NA
        	  }}
          x<-x[names(x)!="Full" & names(x)!=paste("RUN",NbRunEval+1,sep="")] 
          x<-round(mean(x,na.rm=T),4)        
        }

        if(weighting.score =='SomersD'){
          x <- x*2 -1
        }
        }
        
        names(x) <- paste(names(x),y@sp.name,sep=".")
        return(x)}
        ), recursive=T)

  failed <- "none"
  if(sum(is.na(weights))>0){
      
      warning(cat(paste('\n\n\n################################\n',
                        paste('The following bivariate model(s) failed and are not included for building the ESMs\n',
                              sep=" "))),print(names(which(is.na(weights))))) 
      failed <- names(which(is.na(weights)))  
    }

  if(is.null(threshold)){
    if(weighting.score =='AUC'){
      threshold <- 0.5
    }else{
      threshold <- 0
    }}    
    
  if(sum(weights<=threshold &!is.na(weights))>0){
    warning(cat(paste('\n\n\n################################\n',
    paste('The following bivariate model(s) is (are) not included for building the ESMs\n because evaluation score is smaller or equal to the given threshold of',weighting.score,'=',threshold,':\n',
    sep=" "))),print(names(weights[weights <= threshold & !is.na(weights)]))) 
  }
  weights[(weights <= threshold) | is.na(weights)] <- 0 



### BUILD ESMs FOR TESTING DATA
  failed.mod<-"none"
  failed.mod<-lapply(mymodel, function(x){
   # if(x@models.failed!="none"){
     return(if(x@models.failed[1]=="none"){paste(x@sp.name,"none",sep=": ")
     }else{
      x@models.failed
     })
  #}
  })

  test.pred<-lapply(mymodel, function(x){
    x <-get_predictions(x,as.data.frame=T) 
   # x<- x[,grep(paste("RUN",NbRunEval+1,sep=""),colnames(x),invert=T)]
    return(x)
    })

  test.ESM <- NULL
  biva.st2<-do.call(cbind, test.pred)

  for(i in 1:length(models)){
      for(run in 1:(NbRunEval+1)){
        if(length(models)>1){
          test.ESM1 <- apply(biva.st2[, grep(paste("RUN",run,"_",models[i],sep=""),colnames(biva.st2))],1,function(x)weighted.mean(x, weights[grep(models[i],names(weights))],na.rm=T))          
        }else{
          test.ESM1 <- apply(biva.st2[, grep(paste("RUN",run,"_",models[i],sep=""),colnames(biva.st2))],1,function(x)weighted.mean(x, weights,na.rm=T))
        }
        test.ESM <- cbind(test.ESM,test.ESM1)
        colnames(test.ESM)[ncol(test.ESM)] <- paste("RUN",run,"_",models[i],sep="")
      }  
  }
  test.ESM<-as.data.frame(test.ESM)
  DATA<-cbind(1:length(data@data.species), resp.var=data@data.species, test.ESM/1000)  
  if("PA" %in% slotNames(data)){ 
    DATA$resp.var[is.na(DATA$resp.var)] <- 0
    }

### EVALUATE ESMs

EVAL <- NULL
for(i in 1:NbRunEval){
  DATA1 <- cbind(DATA[,1:2],DATA[,grep(paste("RUN",i,"_",sep=""),colnames(DATA))])
  colnames(DATA1)[3] <- colnames(DATA)[2+i]
  if(length(models)>1){
  if(nrow(DATA1) %in%  apply(DATA1[,3:ncol(DATA1)],2,function(x){sum(is.na(x))})){
	failed<- which(nrow(DATA1) ==  apply(DATA1,2,function(x){sum(is.na(x))}))
	DATA1[,nrow(DATA1) ==  apply(DATA1,2,function(x){sum(is.na(x))})]	<- 0
	#DATA1<-DATA1[,nrow(DATA1) !=  apply(DATA1,2,function(x){sum(is.na(x))})]	
  }
  }else{if(nrow(DATA1)==sum(is.na(DATA1[,3]))){
    failed<-  paste("RUN", i, sep = "")
    DATA1[,3] <-0
  }
  }  
  EVAL1<-presence.absence.accuracy(DATA1[!calib.lines[,i],], threshold=as.vector(optimal.thresholds(DATA1[!calib.lines[,i],], opt.methods="MaxSens+Spec"),mode='numeric')[-1])
  EVAL1<-EVAL1[c(1,2,4:7,9:12)]
  EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity-1
  if(length(models)>1){
	if(nrow(DATA1) %in%  apply(DATA1[,3:ncol(DATA1)],2,function(x){sum(is.na(x))})){
  		EVAL1[EVAL1$model==names(failed),2:11]<-NA
	}
  }else{if(nrow(DATA1)==sum(is.na(DATA1[,3]))){
    EVAL1[EVAL1$model==names(failed),2:11]<-NA
  }}
    
  EVAL1$SomersD <- EVAL1$AUC*2-1
  EVAL1$Boyce <- EVAL1$MPA <- NA
  for(n in 1:nrow(EVAL1)){
    EVAL1$MPA[EVAL1$model==EVAL1$model[n]] <- ecospat.mpa(DATA1[!calib.lines[,i] & DATA1[,2] ==1,EVAL1$model[n]])
    # if(packageVersion("ecospat")==1.0){
    #  EVAL1$Boyce[EVAL1$model==EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[,i],EVAL1$model[n]],DATA1[!calib.lines[,i]&DATA1[,2] ==1,EVAL1$model[n]],PEplot = F)$Pearson.cor
    # }else{
      EVAL1$Boyce[EVAL1$model==EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[,i],EVAL1$model[n]],DATA1[!calib.lines[,i]&DATA1[,2] ==1,EVAL1$model[n]],PEplot = F)$Spearman.cor
    #}
  }
  EVAL1$technique <- unlist(strsplit(EVAL1$model,split="_"))[seq(2, nrow(EVAL1)*2, 2)]
  EVAL1$RUN <- paste("RUN",i,sep="")
  EVAL <- rbind(EVAL, EVAL1)
}


### EVALUATE 'DOUBLE ENSEMBLE'
 

if(length(models)>1){ #there is no double-ensemble if only one technique is applied
  
  weights.double <- aggregate(EVAL[,weighting.score],by=list(EVAL$technique), FUN=mean, na.rm=T)
  weights.double <- weights.double[order(weights.double[,1]),]

    if(!is.null(threshold)){
      weights.double[weights.double <= threshold] <- 0
      if(sum(weights.double==0)>0){
        cat(paste('\n\n\n################################\n',
        paste('The following ESM model(s) is (are) not included for building the \'double ensemble\'\n because evaluation score is smaller or equal to the given threshold:\n',
        list(weights.double[weights.double[,2]==0,1]),sep=" ",
        '\n################################\n')))
      }
    }


  for(run in 1:(NbRunEval+1)){
    test.ESM$EF <- apply(test.ESM[,grep(paste("RUN",run,"_",sep=""), colnames(test.ESM))],1,function(x)weighted.mean(x, weights.double[,2],na.rm=T))
    colnames(test.ESM)[ncol(test.ESM)] <- paste("RUN",run,"_EF",sep="")
  }
  DATA<-cbind(1:length(data@data.species), resp.var=data@data.species, test.ESM/1000)  
  if("PA" %in% slotNames(data)){ 
    DATA$resp.var[is.na(DATA$resp.var)] <- 0
  }
  
  EVAL <- NULL
  for(i in 1:NbRunEval){
    DATA1 <- cbind(DATA[,1:2],DATA[,grep(paste("RUN",i,"_",sep=""),colnames(DATA))])
 	 if(nrow(DATA1) %in%  apply(DATA1[,3:ncol(DATA1)],2,function(x){sum(is.na(x))})){
		failed<- which(nrow(DATA1) ==  apply(DATA1,2,function(x){sum(is.na(x))}))
		DATA1[,nrow(DATA1) ==  apply(DATA1,2,function(x){sum(is.na(x))})]	<- 0
		#DATA1<-DATA1[,nrow(DATA1) !=  apply(DATA1,2,function(x){sum(is.na(x))})]	
	  }

    EVAL1<-presence.absence.accuracy(DATA1[!calib.lines[,i],], threshold=as.vector(optimal.thresholds(DATA1[!calib.lines[,i],], opt.methods="MaxSens+Spec"),mode='numeric')[-1])
    EVAL1<-EVAL1[c(1,2,4:7,9:12)]
    EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity-1
	if(nrow(DATA1) %in%  apply(DATA1[,3:ncol(DATA1)],2,function(x){sum(is.na(x))})){
  		EVAL1[EVAL1$model==names(failed),2:11]<-NA
	}
    EVAL1$SomersD <- EVAL1$AUC*2-1
    EVAL1$Boyce <- EVAL1$MPA <- NA
    for(n in 1:nrow(EVAL1)){
      EVAL1$MPA[EVAL1$model==EVAL1$model[n]] <- ecospat.mpa(DATA1[!calib.lines[,i] & DATA1[,2] ==1,EVAL1$model[n]])
      # if(packageVersion("ecospat")==1.0){
      #  EVAL1$Boyce[EVAL1$model==EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[,i],EVAL1$model[n]],DATA1[!calib.lines[,i]&DATA1[,2] ==1,EVAL1$model[n]],PEplot = F)$Pearson.cor
      # }else{
        EVAL1$Boyce[EVAL1$model==EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[,i],EVAL1$model[n]],DATA1[!calib.lines[,i]&DATA1[,2] ==1,EVAL1$model[n]],PEplot = F)$Spearman.cor
      #}
    }
    EVAL1$technique <- unlist(strsplit(EVAL1$model,split="_"))[seq(2, nrow(EVAL1)*2, 2)]
    EVAL1$RUN <- paste("RUN",i,sep="")
    EVAL <- rbind(EVAL, EVAL1)
  }

### BUILD THE 'DOUBLE ENSEMBLE' PROJECTION  

}
  
colnames(DATA) <- gsub(paste("RUN",NbRunEval+1,sep=""),"Full", colnames(DATA))
colnames(DATA)[3:ncol(DATA)] <- paste(colnames(DATA)[3:ncol(DATA)],"ESM",sep="_")
DATA[,-c(1,2)] <- DATA[,-c(1,2)]*1000
EVAL[,2:14] <- round(EVAL[,2:14],3)
EVAL[,2] <- EVAL[,2]*1000

####  FINAL OUTPUT #########################################################

  if(length(models)>1){
  output<-list(species = data@sp.name, ESM.fit=round(DATA[,-1]), ESM.evaluations = EVAL,weights=weights, weights.EF=weights.double, failed= failed.mod)
  }else{
  output<-list(species = data@sp.name, ESM.fit=round(DATA[,-1]), ESM.evaluations = EVAL,weights=weights, failed= failed.mod)
  }
  setwd(iniwd)
   return(output)
}    


############################
# c)  ecospat.ESM.EnsembleProjection    
############################


#Description
## Projecting calibrated ESMs into new space or time

## FUNCTION'S ARGUMENTS
## ESM.prediction.output,:          Object returned by ecospat.ESM.Projection
## ESM.EnsembleModeling.output:     Object returned by ecospat.ESM.EnsembleModeling

## Values:
# ESM.projections:    Returns the projections of ESMs for the selected single models and their ensemble


## Authors:
# Frank Breiner

## References
#Lomba, A., L. Pellissier, C. F. Randin, J. Vicente, F. Moreira, J. Honrado, and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation 143:2647-2657.
# Breiner F.T., Guisan A., Bergamini A., Nobis M.P. (2015):
#  Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6: 1210-1218. doi: 10.1111/2041-210X.12403

## See also
# ecospat.ESM.Modeling; ecospat.ESM.MergeModels; ecospat.ESM.EnsembleModeling

ecospat.ESM.EnsembleProjection <- function(ESM.prediction.output, ESM.EnsembleModeling.output){
  
  models <- ESM.prediction.output$models
  weights <- ESM.EnsembleModeling.output$weights
  weights.EF <- ESM.EnsembleModeling.output$weights.EF
  NbRunEval <- ESM.prediction.output$NbRunEval
  pred.biva <-  ESM.prediction.output$pred.biva
  new.env.raster  <- ESM.prediction.output$new.env.raster
  failed.mod <- grep(paste("RUN",NbRunEval+1,sep=""),unlist(ESM.EnsembleModeling.output$failed),value=T)

  ## to remove weigths for full models which failed!
  
  
  ## MAKE ESM PROJECTIONS
  if(!exists("new.env")){stop("new.env object required!")}
  if(new.env.raster)pred.biva<-grep("gri", pred.biva, value=T)
  if(!new.env.raster)pred.biva<-grep("RData", pred.biva, value=T)
  
  biva.proj<-list()
  for(i in 1:length(pred.biva))
  {
    if(!new.env.raster){
      biva.proj[[i]]<-as.data.frame(get(load(pred.biva[i])))
      colnames(biva.proj[[i]]) <- gsub("AllData", "BIOMOD",colnames(biva.proj[[i]]))
      colnames(biva.proj[[i]]) <- gsub( paste("RUN",NbRunEval+1,sep=""), "ESM",colnames(biva.proj[[i]]))
      colnames(biva.proj[[i]]) <- paste(colnames(biva.proj[[i]]),i, sep=".")
    }
    if(new.env.raster){biva.proj[[i]]<-stack(pred.biva[i])}
  }
  
  ######## Use bivariate models to build ESMs for each model technique 
  
  if(!new.env.raster)
  { 
    pred.ESM <- list()
    biva.st<-do.call(cbind, biva.proj)
    
    for(i in 1:length(models)){
      if(length(models)>1){
        
        ## to remove weigths for full models which failed!
        wm <- weights[grep(models[i],names(weights))][names(weights[grep(models[i],names(weights))]) %in% colnames( biva.st[,grep(models[i],colnames(biva.st))])]
        
        
        pred.ESM[[i]] <- apply(biva.st[,grep(models[i],colnames(biva.st))],1,function(x) stats::weighted.mean(x,wm,na.rm=T))
      }else{
        pred.ESM[[i]] <- apply(biva.st[,grep(models[i],colnames(biva.st))],1,function(x) stats::weighted.mean(x,wm,na.rm=T))      
      }
    };rm(biva.st)
    pred.ESM <- as.data.frame(do.call(cbind, pred.ESM))
    colnames(pred.ESM) <- models
  }
  
  if(new.env.raster){
    pred.ESM<-stack(biva.proj)
    for(i in 1:length(models)){
      assign(models[i], pred.ESM[[grep(models[i], names(pred.ESM))]])
      n <- grep(models[i],failed.mod,value=T)
      if(length(models)>1){
        weights.mod<-weights[grep(models[i], names(weights))]        
        ## Remove weights from models where Full model failed
        if(length(n)>0){
#          weights.mod<- weights.mod[grep( #code by Frank, but doesn't work if more than 1 model failed
#            paste(unlist(strsplit(n,"_"))[((1:length(n))*4)-3],"_",sep=""),
#            paste(names(weights.mod),"_",sep=""), invert=T)]
	   
	   weights.mod<-weights.mod[!names(weights.mod) %in% paste(models[i],".",unlist(strsplit(n,"_"))[((1:length(n))*4)-3],sep="")]	
           assign(models[i], round(raster::weighted.mean(get(models[i]), weights.mod,na.rm=T)) )
        }else{
          assign(models[i], round(raster::weighted.mean(get(models[i]), weights.mod,na.rm=T)) )
        }
      }else{
        if(length(n)>0){
          weights.mod<- weights[grep(
            paste(unlist(strsplit(n,"_"))[1],"_",sep=""),
            paste(names(weights),"_",sep=""), invert=T)]       
          assign(models[i], round(raster::weighted.mean(get(models[i]), weights.mod,na.rm=T)) )
        }else{
          assign(models[i], round(raster::weighted.mean(get(models[i]), weights, na.rm=T)) )        
        }
      }
      
    }
    pred.ESM <- stack(mget(models))[[order(models)]];do.call("rm",as.list(models))
  }
  
  ## Do projection for Double Ensemble
  
  if(length(models)>1){
  if(new.env.raster)
  {    
    ESM.EF  <- round(raster::weighted.mean(pred.ESM[[names(pred.ESM)]], weights.EF[order(weights.EF[,1]),2],na.rm=T))
    pred.ESM <- stack(pred.ESM, ESM.EF)
    names(pred.ESM) <- c(names(pred.ESM)[1:(nlayers(pred.ESM)-1)],"EF") ; rm(ESM.EF)
  }
  
  if(!new.env.raster)
  { 
    pred.ESM$EF <- apply(pred.ESM,1,function(x) stats::weighted.mean(x,weights.EF[order(weights.EF[,1]),2],na.rm=T))
  }
  
  
  if(!new.env.raster)
  {pred.ESM <- round(pred.ESM)}
  }
  
  return(ESM.projections=pred.ESM)
}

############################
# d)  ecospat.mpa   
############################

## This function calculates the Minimal Predicted Area.

## FUNCTION'S ARGUMENTS
## Pred:      numeric or RasterLayer .predicted suitabilities from a SDM prediction
## Sp.occ.xy: xy-coordinates of the species (if Pred is a RasterLayer)
## perc:      Percentage of Sp.occ.xy that should be encompassed by the binary map.

## Details:

## Value
# Returns the Minimal Predicted Area

## Author(s)
# Frank Breiner

## References
# Engler, R., A. Guisan, and L. Rechsteiner. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. Journal of Applied Ecology 41:263-274.


ecospat.mpa <- function(Pred, Sp.occ.xy, perc = 0.9){
  perc <- 1-perc
  if(class(Pred)=="RasterLayer"){
    Pred <- extract(Pred, Sp.occ.xy)}
  round(quantile(Pred,probs=perc),3)}

### EXAMPLE

# obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
#         [which(ecospat.testData$Saxifraga_oppositifolia==1)])
# 
# ecospat.mpa(obs)
# ecospat.mpa(obs,perc=1) ## 100% of the presences encompassed


## Example script for using Ensemble Small Models ESMs according to Lomba et al. 2010
## Written by Frank Breiner
## Swiss Federal Research Institute WSL, July 2014.

## References:
## Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J. & Guisan, A. (2010). Overcoming the rare species modeling paradox: A novel hierarchical framework applied to an Iberian endemic plant. Biological Conservation, 143, 2647-2657.



############################
# e)  ecospat.ESM.MergeModels   merge different ecospat.ESM.Modeling outputs
############################


#Description
# Enables to merge different ecospat.ESM.Modeling outputs which were produced e.g. on different computers or on a cluster


## FUNCTION'S ARGUMENTS
## ESM.modeling.output: list object of different ESM.modeling.output objects. 

#Details
#The Biomod output folders must be copied into one working directory if the output is at different locations!


## Values:
# modeling.id:    "character", id of modeling procedure
# models.:        location of models output
# models:         models argument
# ESM_Projection: logical. TRUE if models were projected for new.env
# pred.biva:      location of the bivariate models predictions
# calib.lines:    "BIOMOD.stored.array", calibrations and evaluation lines selected for each run
# NbRunEval:      Number of Evaluation run
# new.env:        set of explanatory variables onto which models were projected 
# data:           'BIOMOD.formated.data' object
# wd:             working directory
# which.biva:     Which bivariate combinations were used for modeling

## Note: You should use this function only if you want to split your bivariate model runs, e.g. on a clusters.

## See also:
# ecospat.ESM.Modeling,	ecospat.ESM.EnsembleModeling

## Authors:
# Frank Breiner


ecospat.ESM.MergeModels <- function(ESM.modeling.output) 
{
  number.outputs <- length(ESM.modeling.output)
  
  ## tests
  if(length(unique(sapply( ESM.modeling.output,function(x) x[3])))!=1){stop("models of ESM.modeling.output differ")}
  if(length(unique(sapply( ESM.modeling.output,function(x) x[7])))!=1){stop("NbRunEval of ESM.modeling.output differ")}
  if(length(unique(sapply( ESM.modeling.output,function(x) x[8])))!=1){stop("new.env of ESM.modeling.output differ")}
  if(length(unique(sapply( ESM.modeling.output,function(x) x[9])))!=1){stop("data of ESM.modeling.output differ")}          
  #if(length(unique(sapply( ESM.modeling.output,function(x) x[2])))==1){stop("The ESM models you want to merge are all identical")}
  combinations<-combn(colnames(ESM.modeling.output[[1]]$data@data.env.var), 2)
  test<-unique(do.call(c, sapply( ESM.modeling.output,function(x) x[11])),MARGIN=2)
  #miss <- !duplicated(cbind(test,combinations),MARGIN=2)[(ncol(test)+1):(ncol(test)+ncol(combinations))]
  if(ncol(combinations)!=length(test)){
    warning(cat(paste('\n\n\n################################\n',"There are",ncol(combinations)-length(test), "Biva models missing for the ESM. You could run ecospat_ESM_Modeling() with which.biva = ",list(which(!1:ncol(combinations) %in% test)), " to include them.",
                      '\n################################\n',sep=" ")) )
    warning("missing: ",print(combinations[,which(!1:ncol(combinations) %in% test)]))        
  }
  ESM.modeling.output[[1]]$modeling.id <- do.call(c,unique(sapply( ESM.modeling.output,function(x) x[1])))
  ESM.modeling.output[[1]]$models. <- unique(do.call(c,unique(sapply( ESM.modeling.output,function(x) x[2]))))
  ESM.modeling.output[[1]]$pred.biva <- unique(do.call(c,unique(sapply( ESM.modeling.output,function(x) x[5]))))
  ESM.modeling.output[[1]]$which.biva <-  unique(do.call(c,unique(sapply( ESM.modeling.output,function(x) x[11]))))
  
  del.pred <- select.pred <- del.models <- select.models <- NULL 
  for(i in 1:length(ESM.modeling.output[[1]]$which.biva)){
    select.models <- c(select.models,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i],"." , sep=""),   ESM.modeling.output[[1]]$models. ,fixed=T)[1])
    del.models <- c(del.models,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i],"." ,sep=""),   ESM.modeling.output[[1]]$models. ,fixed=T)[-1])
    
    if(length(grep(".grd",   ESM.modeling.output[[1]]$pred.biva ))>0){
      select.pred <- c(select.pred,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i], ".grd",sep=""),   ESM.modeling.output[[1]]$pred.biva )[1])
      select.pred <- c(select.pred,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i], ".gri",sep=""),   ESM.modeling.output[[1]]$pred.biva )[1])
      del.pred <- c(del.pred,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i], ".grd",sep=""),   ESM.modeling.output[[1]]$pred.biva )[-1])
      del.pred <- c(del.pred,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i], ".gri",sep=""),   ESM.modeling.output[[1]]$pred.biva )[-1])     
    }
    if(length(grep(".RData",   ESM.modeling.output[[1]]$pred.biva ))>0){
      select.pred <- c(select.pred,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i], ".RData",sep=""),   ESM.modeling.output[[1]]$pred.biva )[1])
      del.pred <- c(del.pred,grep(paste("ESM.BIOMOD.",  ESM.modeling.output[[1]]$which.biva[i], ".RData",sep=""),   ESM.modeling.output[[1]]$pred.biva )[-1])
    }
  }
  if(length(del.models)>0){ 
    warning(cat(paste('\n\n\n################################\n',
                      length(del.models), "models are not considered for ESMs because they seem to be either replicates or artefacts (e.g. BIOMOD.models.out-objects which already existed before running ecospat_ESM_Modeling function but do not belong to the ESM).\n
                      The following models were deleted:\n",
                      sep=" ")),print(ESM.modeling.output[[1]]$models.[del.models]) ) 
  }
  if(length(del.pred)>0){ 
    warning(cat(paste('\n\n\n################################\n',length(del.pred), "predictions were deleted because they seem to be replicates:\n",
                      sep=" ")),print(ESM.modeling.output[[1]]$pred.biva[del.pred]) )    
  }
  ESM.modeling.output[[1]]$models. <-   ESM.modeling.output[[1]]$models.[select.models]
  ESM.modeling.output[[1]]$pred.biva <-   ESM.modeling.output[[1]]$pred.biva[select.pred]
  
  return(ESM.modeling.output[[1]])
}                                
