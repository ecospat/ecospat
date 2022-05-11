ecospat.ESM.EnsembleEvaluation <- function(ESM.modeling.output,ESM.EnsembleModeling.output,metrics = c("SomersD","AUC","MaxTSS","MaxKappa","Boyce"), EachSmallModels = FALSE){


  #require(dismo) #evaluate function
  
  metrics<-tolower(metrics)
  if(length(intersect(metrics, c("auc", "maxtss", "boyce", "maxkappa", "somersd")))==0){
    stop("The chosen metrics are not supported! Choose at least one of the following: AUC, MaxTSS, Boyce, MaxKappa or SomersD")
  }
  if(!is.logical(EachSmallModels)){
    stop("EachSmallModels should be logical")
  }
  
  resp <- ESM.EnsembleModeling.output$ESM.fit[,1]
  modelling.techniques <- ESM.modeling.output$models
  nMod <- length(ESM.modeling.output$mymodels)
  nReplicate <- ESM.modeling.output$NbRunEval
  fit <- ESM.EnsembleModeling.output$ESM.fit
  calib <- ESM.modeling.output$calib.lines[,1:nReplicate]
  wd <- ESM.modeling.output[["wd"]]
  
  ######################
  
  PredFin <- NULL 
  evalFin <- NULL
  
  for(d in 1:length(modelling.techniques)){
    
    fitMod <- fit[,c(1,grep(modelling.techniques[d],colnames(fit)))] #Select the column for the modelling technique
    fitMod <- fitMod[,-c(grep("Full",colnames(fitMod)))] # Remove the full model
    
    Pred <- .ecospat.pooling(calib=calib,models.prediction=fitMod) 
    
    
    if(d==1){
      PredFin <- cbind(PredFin,Pred)
    }else{
      PredFin <- cbind(PredFin,Pred[,-1])
    }
    
    colnames(PredFin)[ncol(PredFin)] = paste0("Fit_",modelling.techniques[d])
    
    evalInter <- .ecospat.evaluationScores(Pred = Pred,metrics = metrics)
    
    evalFin <- rbind(evalFin,evalInter)
    rownames(evalFin)[nrow(evalFin)] = modelling.techniques[d]
  }
  
  if(length(modelling.techniques)>1){ #If there is more than 1 modelling algorithm, need to evaluate the consensus called here "ensemble"
    
    weights <- ESM.EnsembleModeling.output$weights.EF[,2]
    PredEns <- cbind.data.frame(resp= PredFin[,1],apply(PredFin[,-1], 1, weighted.mean,w=weights))
    PredFin <- cbind(PredFin,PredEns[,-1])
    colnames(PredFin)[ncol(PredFin)] = "Fit_ensemble"
    
    
    ###Computation of the evaluation metrics based on this big data set 
    evalInter <- .ecospat.evaluationScores(Pred = PredEns,metrics = metrics)
    
    evalFin <- rbind(evalFin,evalInter)
    rownames(evalFin)[nrow(evalFin)] = "ensemble"
  }

  output <- list(ESM.evaluations = evalFin, ESM.fit = PredFin)
  
  if(EachSmallModels){ ## Slow take few minutes
    
    evalBivaFin <- list()
    PredBivaFin <- list()
    
    for(i in 1:nMod){
      if(i== round(1*(nMod)/4)){
      }
      if(i== round((nMod)/2)){
      }
      if(i== round(3*(nMod)/4)){
      }
    
      evalBiva <- NULL
      PredBiva <- NULL
      
      IndivMod <- ESM.modeling.output$mymodels[[i]]
      models.prediction <- get(load(paste0(wd,"/",IndivMod@models.prediction@link)))
      models.prediction <- models.prediction[,modelling.techniques[d],1:nReplicate,]
      models.prediction <- cbind.data.frame(resp=resp,models.prediction)
      
      for(d in 1:length(modelling.techniques)){

        Pred <- .ecospat.pooling(calib=calib,models.prediction=models.prediction) 
        PredBiva <- cbind(PredBiva,Pred)
        Pred <- na.omit(Pred) #Remove points where the models failed
        colnames(PredBiva)[ncol(PredBiva)] = paste0("Fit_",modelling.techniques[d])
        
        evalInter <- .ecospat.evaluationScores(Pred = Pred,metrics = metrics)
        
        evalBiva <- rbind(evalBiva,evalInter)
        rownames(evalBiva)[nrow(evalBiva)] = modelling.techniques[d]
      }
      
      evalBivaFin[[i]] = evalBiva
      PredBivaFin[[i]] = PredBiva
    }
    
    output$ESM.evaluations.bivariate.models = evalBivaFin
    output$ESM.fit.bivariate.models = PredBivaFin
  } ## End of each bivariate model evaluations

  return(output)
  
}

.ecospat.pooling <- function(calib,models.prediction){
  Pred <- NULL
  for(k in 1:nrow(calib)){ 
    if(sum(!calib[k,])!=0){ #If the point is used to evaluate at the least one replicate
      valStock <- cbind(models.prediction[k,1],mean(as.numeric(models.prediction[k,(which(!calib[k,]))]),na.rm=T))#if a point is used twice for the evaluation, we take the mean of its fitted values trough the different runs
      colnames(valStock) = c("resp","meanESM")
      Pred <- rbind(Pred,valStock)
    }
  }
  
  return(Pred)
}

.ecospat.evaluationScores <- function(Pred,metrics){
  evalInter <- NULL
  pred.esmPres <-Pred[Pred[,"resp"]==1,2]/1000 #/1000  to have probabilities value
  pred.esmAbs <-Pred[Pred[,"resp"]==0,2]/1000
  
  if("auc" %in% metrics){
    auc.test <- dismo::evaluate(p=pred.esmPres, a=pred.esmAbs)@auc
    evalInter<- cbind(evalInter,AUC=auc.test)
  }
  if("somersd" %in% metrics){
    if("auc" %in% metrics){
      evalInter<- cbind(evalInter,SomersD=(2*auc.test-1))
    }
    else{
      auc.test <- dismo::evaluate(p=pred.esmPres, a=pred.esmAbs)@auc
      evalInter<- cbind(evalInter,SomersD=(2*auc.test-1))
    }
    
  }
  if("boyce" %in% metrics){
    boyce.test <- ecospat.boyce(c(pred.esmPres,pred.esmAbs),pred.esmPres, PEplot=F)$cor
    evalInter <- cbind(evalInter,Boyce=boyce.test)
  }
  if("maxtss" %in% metrics){
    tss.test <- ecospat.max.tss(Pred = (Pred[,2]/1000),Sp.occ = Pred[,1])[[2]]
    evalInter <- cbind(evalInter,MaxTSS=tss.test)
  }
  if("maxkappa" %in% metrics){
    max.kappa.test <- ecospat.max.kappa(Pred = (Pred[,2]/1000),Sp.occ = Pred[,1])[[2]]
    evalInter <- cbind(evalInter,MaxKappa=max.kappa.test)
  }
  return(evalInter)
}


  