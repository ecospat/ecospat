ecospat.nichePOSNB<-function(df,colvar,colfreq){
  
  #require(Hmisc)
  
  if (any(!is.numeric(colvar),colvar==0,is.null(colvar),!colvar%in%1:ncol(df))) stop("colvar should point to relevant columns in df corresponding to environmental axes")
  if (any(!is.numeric(colfreq),colfreq==0,is.null(colfreq),!colfreq%in%1:ncol(df))) stop("colfreq should point to relevant columns in df corresponding to taxa frequencies")
  
  var<-data.frame(df[,colvar])
  freq<-data.frame(df[,colfreq])
  
  POSNB<-matrix(nrow=ncol(freq),ncol=(length(colvar)*2))
  row.names(POSNB)<-colnames(freq)
  colnames(POSNB)<- c(paste0(rep(colnames(df)[colvar],2),rep(c("_pos","_nb"),each=length(colvar))))
  
  for (i in 1:ncol(freq)){
    freq.i<-freq[,i]
    freq.i<-(freq.i-min(freq.i))/(max(freq.i)-min(freq.i))
    rangevar<-1:length(colvar)
    for (j in rangevar){
      POSNB[i,j]<- Hmisc::wtd.mean(var[,j],freq.i)
      POSNB[i,j+length(colvar)] <- sqrt(Hmisc::wtd.var(var[,j],freq.i))
    }
  }
  return(POSNB)
}

ecospat.nicheNBmean<-function(POSNB,w=NULL){
  if(ncol(POSNB)==2) stop("POSNB should have more than one axis")
  NB<-POSNB[,(ncol(POSNB)/2+1):ncol(POSNB)]
  if(is.null(w)) w<-rep(1/ncol(NB),ncol(NB))
  if (!is.null(w)&length(w)!=ncol(NB)) stop("w should have the same length as the number of axes in POSNB")
  NBmean<-NB%*%w/sum(w)
  colnames(NBmean)<-"NBmean"
  return(NBmean)
  }


## test
#df <- read.delim("ecospat/data/ecospat.testNichePOSNB.txt")
#POSNB<-ecospat.nichePOSNB(df,colvar=c(2),colfreq = 4:15) # 1 axes
#POSNB<-ecospat.nichePOSNB(df,colvar=c(2:3),colfreq = 4:15) # 2 axes
#POSNB<-ecospat.nichePOSNB(df,colvar=c(2:5),colfreq = 4:15) # 4 axes
#ecospat.nicheNBmean(POSNB,w=c(2,2))

