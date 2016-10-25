#######################################################################################################
## FUNCTION TO APPLY IN THE SESAM FRAMEWORK TO MODEL COMMUNITY COMPOSITION
  
##Probability Ranking Rule 

##M. D`Amen, Dept.of Ecology and Evolution, University of Lausanne    
######################################################################################################
#
# Function to implement the SESAM framework with the `Probability Ranking` Rule as biotic rule.
# The community composition in each site is determined by ranking the species in decreasing order of their 
# predicted probability of presence from SDMs (other criteria to assign probabilities to the species 
# can be applied) up to a richness prediction.
# For further explanations see:  
# D`Amen et al., 2015. Using species richness and functional traits predictions to constrain assemblage 
# predictions from stacked species distribution models. J Biogeogr. doi:10.1111/jbi.12485

# proba: probabilities from SDMs (or other sources) for all species (need to have defined row names)
# sr: richness predictions (data frame with richness value in the first column) 
# the function saves in the working directory a .txt file with the community prediction by the SESAM framework

################################
####Probability Ranking Rule ###
################################

ecospat.SESAM.prr <- function(proba, sr) {
  
  projSR <- round(round(as.vector(sr[[1]])))
  
  new.prob.prr <- proba
  dataSSDM_p <- proba
  
  for (i in 1:nrow(proba)) {
    
    print(paste("test.prr, processing row ", i, sep = ""))
    
    SR <- projSR[i]  # values of species richeness     
    if (SR > 0) {
      predcom <- dataSSDM_p[i, ]
      predcom_p <- dataSSDM_p[i, ]  #corresponding row with HS 
      com <- order(predcom_p, decreasing = TRUE)
      pres <- com[1:SR]
      predcom[, pres] <- 1
      predcom[, -pres] <- 0
    } else {
      predcom[, ] <- 0
    }
    
    new.prob.prr[i, ] <- predcom
  }
  
  # return(community.prediction.prr<-new.prob.prr)
  print(new.prob.prr)
  
  write.table(new.prob.prr, "community.prediction.prr.txt", sep = "\t")  #final outcome of SESAM: community prediction
}