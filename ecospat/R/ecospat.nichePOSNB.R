
ecospat.nichePOSNB<-function(df,colvar,colfreq){
  
  require(Hmisc)
  var<-data.frame(df[,colvar])
  freq<-data.frame(df[,colfreq])
  
  if (!is.numeric(colvar)||colvar==0||is.null(colvar)) stop("colvar should point to relevant columns in df corresponding to environmental axes")
  if (!is.numeric(colfreq)||colfreq==0||is.null(colfreq)) stop("colfreq should point to relevant columns in df corresponding to taxa frequencies")
  
  if(ncol(var)==1) {
    meanNicheByTaxa<-matrix(nrow=ncol(freq),ncol=2)
    row.names(meanNicheByTaxa)<-colnames(freq)
    colnames(var)<-colnames(df)[colvar]
    colnames(meanNicheByTaxa)<- c(paste0(colnames(var),rep(c("_pos","_nb"))))
  }
  
  if(ncol(var)!=1){
    meanNicheByTaxa<-matrix(nrow=ncol(freq),ncol=(length(colvar)*2)+2)
    row.names(meanNicheByTaxa)<-colnames(freq)
    colnames(meanNicheByTaxa)<- c(paste0(rep(colnames(var),2),rep(c("_pos","_nb"),each=length(colvar))),"mean_pos","mean_nb")
  }
  
  for (i in 1:nrow(meanNicheByTaxa)){
    freq.i<-freq[,i]
    freq.i<-(freq.i-min(freq.i))/(max(freq.i)-min(freq.i))
    rangevar<-1:length(colvar)
    for (j in rangevar){
      meanNicheByTaxa[i,j]<- Hmisc::wtd.mean(var[,j],freq.i)
      meanNicheByTaxa[i,j+length(colvar)] <- sqrt(Hmisc::wtd.var(var[,j],freq.i))
    }
    if(ncol(var)!=1){
      meanNicheByTaxa[i,ncol(meanNicheByTaxa)-1] <- mean(meanNicheByTaxa[i,rangevar])
      meanNicheByTaxa[i,ncol(meanNicheByTaxa)] <- mean(meanNicheByTaxa[i,rangevar+length(colvar)])
    }
  }
  return(meanNicheByTaxa)
}

## test
<<<<<<< HEAD
#df <- read.delim("ecospat/data/ecospat.testNichePOSNB.txt")
#ecospat.nichePOSNB(df,colvar=c(2),colfreq = 4:15) # 1 axes
#ecospat.nichePOSNB(df,colvar=c(2:3),colfreq = 4:15) # 2 axes
#ecospat.nichePOSNB(df,colvar=c(2:6),colfreq = 4:15) # 5 axes
=======
# df <- read.delim("ecospat/data/ecospat.testNichePOSNB.txt")
# df <-data(ecospat.testNichePOSNB)
# ecospat.nichePOSNB(df,colvar=c(2),colfreq = 4:15) # 1 axes
# ecospat.nichePOSNB(df,colvar=c(2:3),colfreq = 4:15) # 2 axes
# ecospat.nichePOSNB(df,colvar=c(2:6),colfreq = 4:15) # 5 axes


>>>>>>> 3cbd64cd7c8e64c36a5b66738741d5ba0294f225
