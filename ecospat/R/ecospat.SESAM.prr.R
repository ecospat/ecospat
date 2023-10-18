#######################################################################################################
## FUNCTION TO APPLY IN THE SESAM FRAMEWORK TO MODEL COMMUNITY COMPOSITION

##Probability Ranking Rule 

##Modified from M. D`Amen, Dept.of Ecology and Evolution, University of Lausanne    
# by V. Verdon, Dept.of Ecology and Evolution, University of Lausanne    
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
# sr: (optional) data frame with richness target value in the first column (output from richness models). 
# If not provided, the richness value will correspond to the sum of the predicted probabilities of all species (see D`Amen et al., 2015.)
# In the rare case when several species have the exact same predicted probabilites and the target richness does not allow all of them to be included 
# within the predicted community, selection between these species is done at random.
# the function saves in the working directory a ".txt" file with the community prediction by the SESAM framework

################################
####Probability Ranking Rule ###
################################


ecospat.SESAM.prr <- function(proba,sr=NULL , verbose = FALSE) {
  # Get target richness per site, either use input richness, or calculate it using sum of probabilities
  if (is.null(sr)){
    row_sums <- round(round(as.vector(rowSums(proba))))
  } else{ 
    row_sums <-round(round(as.vector(sr[[1]])))
  }
  
  # Sort the probabilities within each row in decreasing order
  sorted_probs <- t(apply(as.matrix(proba), 1, function(x) sort(x, decreasing = TRUE)))
  
  # Get threshold probability at which max richness of each site is reached
  threshold<-sorted_probs[cbind(1:nrow(sorted_probs),row_sums)]
  
  #binarize probability matrix using site specific threshold
  binary_dtframe <- 1*data.frame(proba>=threshold)

  #correction when several species with same predicted probability have been kept within community exceeding input richness.
  new_rowsums<-rowSums(binary_dtframe)
  diff_rowsums<-new_rowsums-row_sums
  if (any(diff_rowsums!=0)){
    positions_potentially_removed<-1*data.frame(proba==threshold)
    positions_toremove <-   as.data.frame(t(mapply(function(row,diff_rowsums) {
      one_indices <- which(row == 1)
      if (diff_rowsums==0){
        row[] <- 0
      } else {
        random_index <- sample(one_indices, diff_rowsums,replace=FALSE)
        row[-random_index] <- 0
      }
      return(as.numeric(row))
    }
    , split(binary_dtframe, 1:nrow(binary_dtframe)), diff_rowsums)))
    

    binary_dtframe<-binary_dtframe*(!(positions_toremove))
  }
  
  return(binary_dtframe)
  
}

