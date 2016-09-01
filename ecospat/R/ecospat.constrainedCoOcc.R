###########################################################
###########################################################
#########        	CO-OCCURRENCE ANALYSIS			     ########
#########  					         &         					   ########
#########  ENVIRONMENTALLY CONSTRAINED NULL MODELS ########
#########                  CtRA1                   ########
###########################################################
###########################################################
# 
# Adapted by Anne Dubuis from P.R Peres-Neto codes used for: 
#
# Peres-Neto P.R., Olden J.D. and Jackson D.A. (2001) Environmentally constrained null models: sites suitablity as occupancy criterion. 
# Oikos 93:110-120
#
# and modified by Manuela D'Amen.
#
# Constrperm() produces a null matrix of 0/1 in which the 1 are weighted by probability values
# from the db "proba". The sum column is costant i.e. null species have the same occurrences as the observed ones
# while row sum is free: richness is variable in the null community sites.
#
# SpeciesCooccurrenceStats() produces a classical C-score index for any confusion matrix: this is used to calculate
# both the C-score for the observed dataset and for the null matrices from  the former functions
#
# NullModels() uses the former two functions and returns the C-score index for the observed community (ObsCscoreTot), the mean of C-score for the simulated communities (SimCscoreTot)
# p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the path specified where the same 
# metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version)
# # NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002 - Ecology)
#
# presence: presence absence table
# pred: table with values of probability of presence for each species in each site
# nbpermut: number of permutation in the null model
# outpath: path to specify where to save the tables
# NB: Format required for imput databases: a plots (rows) x species (columns) matrix
# Input matrices should have column names (species names) and row names (sampling plots)


##################################################################################################################################################

# NullModels(presence, pred, 1000, outpath) # launch command
# library(ade4)

#############################
## Constrained permutation ##
#############################
Constrperm<-function(presence, pred) 
{
  nbsps<-ncol(presence) # species number
  nbsites<-nrow(presence) # sites number
  
  nbocc<- as.vector(apply(presence,MARGIN=2,sum)) # occurrences number for each species
  sumprob<-as.vector(apply(pred,MARGIN=2,sum)) # occurrence probability sum for each species
  
  matsum<-matrix(0,ncol=nbsps,nrow=nbsites)
  
  for(i in 1:nbsps)
  {
    matsum[,i]<-sumprob[i] 
  }
  
  # creates a relative probability matrix
  pred<-pred/matsum  
  
  transpo<-t(as.matrix(presence))
  sps.names<-row.names(transpo)
  sites<-row.names(presence)
  noms<-list(sites, sps.names)
  
  # creates a random matrix weighted by the species' probabilities of occurrence
  nullmat<-matrix (0,nrow=nbsites,ncol=nbsps,dimnames=noms) 
  randr<-matrix(runif(nbsites*nbsps),nbsites)
  randr<-randr*pred 
  
  for (i in 1: nbsps) # matrix randomisation
  {
    a<-as.vector(randr[,i])
    names(a)<-(1:nbsites)
    a<-sort(a)
    b<-as.numeric(names(a))
    x<-nbsites-nbocc[i]+1
    
    for(j in x:nbsites){
      r<-b[j]
      nullmat[r,i]<-1
    }
  }
  return (nullmat)
}

### END function ###


##############################
## Co-occurrence statistics ##
##    calculations          ##
##############################
SpeciesCooccurrenceStats<-function (presence)  
{
  nbsps<-ncol(presence) # species number
  nbsites<-nrow(presence) # sites number
  nbocc<- as.vector(apply(presence,MARGIN=2,sum))	# occurrences number for each species
  
  # C-score for each species
  presence<-as.matrix(presence)
  coocc<-t(presence)%*%presence
  nbspec=dim(coocc)[1]
  mat1<-array(apply(presence,MARGIN=2,sum),dim=c(nbsps,nbsps))
  mat2<-t(array(apply(presence,MARGIN=2,sum),dim=c(nbsps,nbsps)))
  Cscoreperspecies <- (mat1 - coocc)*(mat2 - coocc)
  
  # C-score on all matrix
  df.synthesis1 <- data.frame(Col = rep(1:ncol(Cscoreperspecies),each=ncol(Cscoreperspecies)),
                              Row = rep(1:nrow(Cscoreperspecies),nrow(Cscoreperspecies)),Sps1 = rep(colnames(Cscoreperspecies),
                                                                                                    each=ncol(Cscoreperspecies)), Sps2 = rep(rownames(Cscoreperspecies),nrow(Cscoreperspecies)),CScore= c(Cscoreperspecies)) 
  v.diago.inf <- c(rownames(df.synthesis1)[df.synthesis1[,1]>df.synthesis1[,2]],rownames(df.synthesis1)[df.synthesis1[,1]==df.synthesis1[,2]])
  df.synthesis <- df.synthesis1[-as.numeric(v.diago.inf),]
  Cscore<-mean(df.synthesis[,5]) 
  
  return(list(coocc=Cscoreperspecies, synth=df.synthesis1, synth_lt=df.synthesis, Cscore=Cscore))
  
}### END function ###


#################
## Null models ##
#################

ecospat.cons_Cscore<- function(presence,pred,nbpermut,outpath)
{ 	
  cat("Computing observed co-occurence matrix", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  
  nbsps<-ncol(presence) 
  nbsites<-nrow(presence) 
  
  # Co-occurrence statistics on observed matrix
  Obs<-SpeciesCooccurrenceStats(presence) 
  
  CooccObs<-Obs$coocc 
  synthObs<-Obs$synth_lt #summary table with Co-occurrence values for each couple
  CscoreTot<-Obs$Cscore
  
  CooccProb<-matrix(0,nrow=nrow(synthObs),ncol=nbpermut, dimnames = list(c(paste(synthObs[,3],synthObs[,4])),c(1:nbpermut)) )
  
  cat("Computing permutations", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  
  if (nbpermut>0){
    for (z in 1:nbpermut){
      #cat(z, "\n",append =TRUE)
      degenerated<-TRUE
      while (degenerated==TRUE){ # looking for a matrix without degenerated sites
        random.distrib.matrix<-Constrperm(presence,pred) 
        sumSites<- as.vector(apply(random.distrib.matrix,MARGIN=2,sum)) 
        for (i in 1:nbsites){
          ifelse(sumSites[i]==0,degenerated<-TRUE,degenerated<-FALSE)	
        }
      }
      
      # Co-occurrence statistics on randomized matrix
      Rnd<-SpeciesCooccurrenceStats(random.distrib.matrix) 
      
      synthRnd<-Rnd$synth_lt	
      
      CooccProb[,z] <- synthRnd[,5] # pairwise C-score for simulated matrix 
    }
    
    ## for the whole community
    vec.CScore.tot<-as.vector(apply(CooccProb,MARGIN=2,mean)) # C-score for all null communities (mean on the columns)
    SimulatedCscore<-mean(vec.CScore.tot) # mean of Simulation C-score: Simulated C-score
    sd.SimulatedCscore<-sd(vec.CScore.tot) # standard deviation of null communities
    ses<-(CscoreTot-SimulatedCscore)/sd.SimulatedCscore # standardized effect size: It scales the results in units of standard deviations, 
    # which allows for meaningful comparisons among different tests
    
    randtest.less<-as.randtest(vec.CScore.tot, CscoreTot, alter="less")
    pval.less<-randtest.less$pvalue
    randtest.greater<-as.randtest(vec.CScore.tot, CscoreTot, alter="greater")
    pval.greater<-randtest.greater$pvalue
    plot(randtest.greater, xlab= "Simulated C-scores",main=paste("", sep=""))
    
    ## for species pairs
    vec.Cscore.pairs<-as.vector(apply(CooccProb,MARGIN=1,mean)) # Mean of null communities C-score for any pair of species
    mat_pval <- matrix(0,nrow(synthRnd),5) 
    mat_pval[,1] <- synthObs[,5] # Observed C-score in the first column
    mat_pval[,2] <- vec.Cscore.pairs # Mean of simulated C-scores in the second column
    
    for (i in 1:nrow(CooccProb))
    {		  
      randtest.less<-as.randtest(CooccProb[i,], mat_pval[i,1], alter="less")
      mat_pval[i,3]<-randtest.less$pvalue
      randtest.greater<-as.randtest(CooccProb[i,], mat_pval[i,1], alter="greater")
      mat_pval[i,4]<-randtest.greater$pvalue
      mat_pval[i,5]<-(mat_pval[i,1]-mean(CooccProb[i,]))/(sd(CooccProb[i,])) #ses for any pair of species
    }    
  }	
  
  cat("Permutations finished",date(), "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat("Exporting dataset", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  cat(".............", "\n",append = FALSE)
  
  df.synthesis <- data.frame(Col = rep(1:ncol(CooccObs),each=ncol(CooccObs)),Row = rep(1:nrow(CooccObs),
                                                                                       nrow(CooccObs)),Sps1 = rep(colnames(CooccObs),each=ncol(CooccObs)), 
                             Sps2 = rep(rownames(CooccObs),nrow(CooccObs)),Cooc_score = c(CooccObs)) 	
  v.diago.inf <- c(rownames(df.synthesis)[df.synthesis[,1]>df.synthesis[,2]],rownames(df.synthesis)[df.synthesis[,1]==df.synthesis[,2]])
  df.synthesis <- df.synthesis[-as.numeric(v.diago.inf),]
  df.synthesis<- cbind(df.synthesis,mat_pval[,2:5])
  names(df.synthesis)[5:9]<-c("C-scoreObs","C-scoreExp","p.less","p.greater","ses")
  
  hist(as.vector(df.synthesis[,9]), xlab="ses", main = paste("Histogram of standardized effect size"))
  abline(v=c(2,-2),col = "red")
  #write.table(df.synthesis,file=paste(outpath,"Results_const_Cscores.txt", sep=""),sep="\t",append= FALSE,row.names= FALSE,col.names=TRUE,quote=FALSE)
  
  #selection of rows with <=0.05
  tab<-df.synthesis
  v<-c(0)
  for (i in 1:nrow(tab)){
    if (tab[i,7]<=0.05||tab[i,8]<=0.05){
      v<-c(v,i)
    }
  }
  m<-data.frame()
  for(j in 1:length(v)){
    m<-rbind(m,tab[v[j],])
  }
  
  m1<-na.omit(m)
  
  write.table(m1,file=paste(outpath,"/Signific_const_Cscores.txt",sep=""),sep="\t",append= FALSE,row.names= FALSE,col.names=TRUE,quote=FALSE)
  
  l<-list(ObsCscoreTot=CscoreTot, SimCscoreTot=SimulatedCscore, PVal.less=pval.less,PVal.greater=pval.greater,SES.Tot=ses)
  return(l)
  
}
### END function ###
