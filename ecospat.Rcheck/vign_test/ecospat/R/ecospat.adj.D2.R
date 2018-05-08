ecospat.adj.D2.glm <- function(glm.obj)
############
## ADJ.D2 ##
############

# Function originally written by A. Guisan (UNIL-ECOSPAT,Switzerland) and modified by C. Randin
# (UNIL-ECOSPAT,Switzerland) S-Plus Function for calculating an adjusted D2 (see Weisberg, 1980) Takes any GLM
# object as argument

{
  go <- glm.obj
  D2 <- (go$null.deviance - go$deviance)/go$null.deviance
  p <- length(go$coefficients)
  n <- length(go$fitted)
  adj.D2 <- 1 - ((n - 1)/(n - p)) * (1 - D2)
  if (adj.D2 < 0) {
    adj.D2 <- 0
    return(adj.D2)
  } else {
    return(adj.D2)
  }
}