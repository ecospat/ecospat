ecospat.poolingEvaluation <- function(fit,calib,resp,AlgoName = NULL,metrics = c("SomersD","AUC","MaxTSS","MaxKappa","Boyce"),ensembleEvaluation=FALSE,w=NULL,metricToEnsemble = "MaxTSS"){
  
  
  #require(dismo) #evaluate function
  
  metrics<-tolower(metrics)
  metricToEnsemble <- tolower(metricToEnsemble)
  if(length(intersect(metrics, c("auc", "maxtss", "boyce", "maxkappa", "somersd")))==0){
    stop("The chosen metrics are not supported! Choose at least one of the following: AUC, MaxTSS, Boyce, MaxKappa or SomersD")
  }
  if(!is.numeric(resp)){
    stop("resp should be numerical vector with 1 for presences and 0 for absences (or background points, pseudo absences) ")
  }
  
  if(!is.logical(as.matrix(calib))){
    stop("calib should be a logical matrix")
  }
  
  if(length(resp)!=nrow(calib) | length(resp)!=nrow(fit[[1]])){
    stop("the length of resp does not match with the number of rows of calib and/or fit")
  }
  if(ncol(calib)!= ncol(fit[[1]])){
    stop("calib and the elements of fit should have the same number of column")
  }
  nAlgo = length(fit)
  if(ensembleEvaluation & !is.null(w) & length(w)!= nAlgo){
    stop("w should have a length corresponding to the value of nAlgo in order to test the ensemble models")
  }
  if(ensembleEvaluation & is.null(w)){
    if(length(metricToEnsemble)>1){stop("metricToEnsemble should have only one element")}
    if(!(metricToEnsemble %in% metrics)){stop("metricToEnsemble should also be in the metrics object")}
    
  }
  if(is.null(AlgoName)){
    AlgoName = 1:nAlgo 
  }
  
  ######################
  
  PredFin <- NULL 
  evalFin <- NULL
  
  for(d in 1:nAlgo){
    
    fitMod <- fit[[d]]
    fitMod <-cbind(resp=resp,fitMod)
    Pred <- .ecospat.pooling(calib=calib,models.prediction=fitMod) 
    
    
    if(d==1){
      PredFin <- cbind(PredFin,Pred)
    }else{
      PredFin <- cbind(PredFin,Pred[,-1])
    }
    
    colnames(PredFin)[ncol(PredFin)] = paste0("Fit_",AlgoName[d])
    
    evalInter <- .ecospat.evaluationScores(Pred = Pred,metrics = metrics)
    
    evalFin <- rbind(evalFin,evalInter)
    rownames(evalFin)[nrow(evalFin)] = AlgoName[d]
  }
  
  if(ensembleEvaluation){
    if(is.null(w)){
      if("auc" %in% metricToEnsemble){metricToEnsemble = "AUC"}
      if("somersd" %in% metricToEnsemble){metricToEnsemble="SomersD"}
      if("boyce" %in% metricToEnsemble){metricToEnsemble = "Boyce"}
      if("maxtss" %in% metricToEnsemble){metricToEnsemble = "MaxTSS"}
      if("maxkappa" %in% metricToEnsemble){metricToEnsemble = "MaxKappa"}
      w <-as.numeric(evalFin[,metricToEnsemble])
    }
    PredEns <- cbind.data.frame(resp = PredFin[,1],apply(PredFin[,-1], 1, weighted.mean,w=w))
    PredFin <- cbind(PredFin,PredEns[,-1])
    colnames(PredFin)[ncol(PredFin)] = "Fit_ensemble"
    
    
    ###Computation of the evaluation metrics based on this big data set 
    evalInter <- .ecospat.evaluationScores(Pred = PredEns,metrics = metrics)
    
    evalFin <- rbind(evalFin,evalInter)
    rownames(evalFin)[nrow(evalFin)] = "ensemble"
  }
  
  output <- list(evaluations = evalFin, fit = PredFin)
  
  return(output)
  
}

ecospat.ESM.EnsembleEvaluation <- function(ESM.modeling.output,ESM.EnsembleModeling.output,metrics = c("SomersD","AUC","MaxTSS","MaxKappa","Boyce"), EachSmallModels = FALSE){
  metrics <- tolower(metrics)
  if (length(intersect(metrics, c("auc", "maxtss", "boyce", 
                                  "maxkappa", "somersd"))) == 0) {
    stop("The chosen metrics are not supported! Choose at least one of the following: AUC, MaxTSS, Boyce, MaxKappa or SomersD")
  }
  if (!is.logical(EachSmallModels)) {
    stop("EachSmallModels should be logical")
  }
  resp <- ESM.EnsembleModeling.output$ESM.fit[, 1]
  modelling.techniques <- ESM.modeling.output$models
  nMod <- length(ESM.modeling.output$mymodels)
  nReplicate <- ESM.modeling.output$NbRunEval
  fit <- ESM.EnsembleModeling.output$ESM.fit
  calib <- ESM.modeling.output$calib.lines[, 1:nReplicate]
  wd <- ESM.modeling.output[["wd"]]
  PredFin <- NULL
  evalFin <- NULL
  for (d in 1:length(modelling.techniques)) {
    fitMod <- fit[, c(1, grep(modelling.techniques[d], colnames(fit)))]
    fitMod <- fitMod[, -c(grep("Full", colnames(fitMod)))]
    Pred <- .ecospat.pooling(calib = calib, models.prediction = fitMod)
    Pred[, -1] = (Pred[, -1])/1000
    if (d == 1) {
      PredFin <- cbind(PredFin, Pred)
    }
    else {
      PredFin <- cbind(PredFin, Pred[, -1])
    }
    colnames(PredFin)[ncol(PredFin)] = paste0("Fit_", modelling.techniques[d])
    evalInter <- .ecospat.evaluationScores(Pred = Pred, metrics = metrics)
    evalFin <- rbind(evalFin, evalInter)
    rownames(evalFin)[nrow(evalFin)] = modelling.techniques[d]
  }
  if (length(modelling.techniques) > 1) {
    weights <- ESM.EnsembleModeling.output$weights.EF[, 2]
    PredEns <- cbind.data.frame(resp = PredFin[, 1], apply(PredFin[, 
                                                                   -1], 1, weighted.mean, w = weights))
    PredFin <- cbind(PredFin, PredEns[, -1])
    colnames(PredFin)[ncol(PredFin)] = "Fit_ensemble"
    evalInter <- .ecospat.evaluationScores(Pred = PredEns, 
                                           metrics = metrics)
    evalFin <- rbind(evalFin, evalInter)
    rownames(evalFin)[nrow(evalFin)] = "ensemble"
  }
  output <- list(ESM.evaluations = evalFin, ESM.fit = PredFin)
  if (EachSmallModels) {
    evalBivaFin <- list()
    PredBivaFin <- list()
    for (i in 1:nMod) {
      evalBiva <- NULL
      PredBiva <- NULL
      IndivMod <- ESM.modeling.output$mymodels[[i]]
      models.File <- get(load(paste0(wd, "/", IndivMod@models.prediction@link)))
      
      for (d in 1:length(modelling.techniques)) {
        models.prediction <- matrix(models.File$pred[models.File$run !="allRun" & models.File$algo==modelling.techniques[d]],ncol=nReplicate,nrow=length(resp),byrow = F)
        models.prediction <- cbind.data.frame(resp = resp, 
                                              models.prediction)
        Pred <- .ecospat.pooling(calib = calib, models.prediction = models.prediction)
        Pred[, -1] = (Pred[, -1]/1000)
        PredBiva <- cbind(PredBiva, Pred)
        Pred <- na.omit(Pred)
        colnames(PredBiva)[ncol(PredBiva)] = paste0("Fit_", 
                                                    modelling.techniques[d])
        evalInter <- .ecospat.evaluationScores(Pred = Pred, 
                                               metrics = metrics)
        evalBiva <- rbind(evalBiva, evalInter)
        rownames(evalBiva)[nrow(evalBiva)] = modelling.techniques[d]
      }
      evalBivaFin[[i]] = evalBiva
      PredBivaFin[[i]] = PredBiva
    }
    output$ESM.evaluations.bivariate.models = evalBivaFin
    output$ESM.fit.bivariate.models = PredBivaFin
  }
  return(output)
}


.ecospat.pooling <- function(calib,models.prediction){
  Pred <- NULL
  for(k in 1:nrow(calib)){ 
    if(sum(!calib[k,])!=0){ #If the point is used to evaluate at the least one replicate
      valStock <- cbind(models.prediction[k,1],mean(as.numeric(models.prediction[k,(which(!calib[k,])+1)]),na.rm=T))#if a point is used twice for the evaluation, we take the mean of its fitted values trough the different runs
      colnames(valStock) = c("resp","meanESM")
      Pred <- rbind(Pred,valStock)
    }
  }
  
  return(Pred)
}

.ecospat.evaluationScores <- function(Pred,metrics){
  evalInter <- NULL
  pred.esmPres <-Pred[Pred[,"resp"]==1,2] 
  pred.esmAbs <-Pred[Pred[,"resp"]==0,2]
  
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
    tss.test <- ecospat.max.tss(Pred = (Pred[,2]),Sp.occ = Pred[,1])[[2]]
    evalInter <- cbind(evalInter,MaxTSS=tss.test)
  }
  if("maxkappa" %in% metrics){
    max.kappa.test <- ecospat.max.kappa(Pred = (Pred[,2]),Sp.occ = Pred[,1])[[2]]
    evalInter <- cbind(evalInter,MaxKappa=max.kappa.test)
  }
  return(evalInter)
}
