ecospat.ESM.threshold <- function(ESM.EnsembleModeling.output, PEplot=FALSE){

  Full.models <- grep('Full',colnames(ESM.EnsembleModeling.output$ESM.fit),value = TRUE)
  
  EVAL <- NULL
  for(i in Full.models){
DATA <- cbind(1:nrow(ESM.EnsembleModeling.output$ESM.fit),
              resp.var = ESM.EnsembleModeling.output$ESM.fit$resp.var, 
              ESM.EnsembleModeling.output$ESM.fit[,i]/1000)


EVAL1 <- presence.absence.accuracy(DATA[, ], threshold = as.vector(optimal.thresholds(DATA[, ], opt.methods = "MaxSens+Spec"), mode = "numeric")[-1])
TSS.th <- EVAL1$threshold
EVAL1 <- EVAL1[c(1, 4:7, 9:12)]
EVAL1$SomersD <- EVAL1$AUC * 2 - 1
boyce <- ecospat.boyce(DATA[,3], DATA[DATA[, 2] == 1,3],PEplot=PEplot)
EVAL1$Boyce <- boyce$Spearman.cor
EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 1
EVAL1$TSS.th <- TSS.th
EVAL1$Boyce.th.max <- EVAL1$Boyce.th.min <- EVAL1$Boyce <- EVAL1$MPA0.90 <-EVAL1$MPA0.95 <-EVAL1$MPA1.0 <- NA

  EVAL1$MPA1.0 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=1)
  EVAL1$MPA0.95 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=.95)
  EVAL1$MPA0.90 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=.9)



  pos.F <- which(boyce$F.ratio>1)
  neg.F <- which(boyce$F.ratio<=1)
  if(max(neg.F) < min(pos.F)){
    EVAL1$Boyce.th.max <- EVAL1$Boyce.th.min <- mean(boyce$HS[c(max(neg.F),min(pos.F))])
    
  }else{
    EVAL1$Boyce.th.max <- mean(boyce$HS[c(max(neg.F),max(neg.F)+1)])
    EVAL1$Boyce.th.min <- mean(boyce$HS[c(min(pos.F),min(pos.F)-1)])
  }
  EVAL1$model <- i
  EVAL <- rbind(EVAL,EVAL1)
  }
  return(EVAL)
}
