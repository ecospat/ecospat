ecospat.co_occurrences <- function (data)
{

#######################
## Co-occurrences 2m ##
#######################

  #---------------------------------------------------------------------------------
  # This part calculate the matrix P1 of co-occurrence realized on co-occurrences
  # possible at a resolution of 2 meters for each pair of species. These results
  # are also put in a vector L1, where each comparison is only kept once.
  #---------------------------------------------------------------------------------
  
  Nsp1 <- ncol(data) - 4
  Nrel1 <- nrow(data)
  P1_n_occ <- t(as.matrix(data[, -(1:4)])) %*% as.matrix(data[, -(1:4)])
  Tot1_sp <- apply(data[, -(1:4)], MARGIN = 2, sum)
  Min1_n <- outer(Tot1_sp, Tot1_sp, FUN = pmin)
  P1_n <- P1_n_occ/Min1_n
  
  #----------------------------------------------------------------------------------
  # This part reunite in columns all the values of co-occurrences relative to each
  # species. Hence, one comparison between two pairs of species occur in two
  # columns.
  #----------------------------------------------------------------------------------
  
  L1_n <- P1_n[lower.tri(P1_n)]
  Sp1_n <- rep(1:(Nsp1 - 1), times = ((Nsp1 - 1):1))
  Sp2_n <- t(outer(1:Nsp1, 1:Nsp1, FUN = "pmax"))[lower.tri(t(outer(1:Nsp1, 1:Nsp1, 
    FUN = "pmax")))]
  # Sumdata_n = data.frame(cbind(Sp1_n,Sp2_n,L1_n))
  
  #----------------------------------------------------------------------------------
  # This part below reunite in columns all the values of co-occurrences relative to
  # each species. There is one column by species. Hence, one comparison between two
  # pairs of species occur in two columns.
  #----------------------------------------------------------------------------------
  
  CoobySp1_n <- matrix(data = 0, nrow = Nsp1 - 1, ncol = Nsp1)
  CoobySp1_n[lower.tri(CoobySp1_n, diag = TRUE)] <- P1_n[lower.tri(P1_n)]
  CoobySp1_n[upper.tri(CoobySp1_n)] <- P1_n[upper.tri(P1_n)]
  CoobySp1_n <- data.frame(cbind(CoobySp1_n))
  # colnames(CoobySp1_n)<-Names
  boxplot(CoobySp1_n)
  
  return(P1_n)
}


## Pairwise CO-OCCURENCE ANALYSIS with C-score  calculation  ##

#C. RANDIN and M. D'Amen, Dept.of Ecology & Evolution, University of Lausanne

## Format required: a plots (rows) x species (columns) matrix of presences/absences
## Input matrix should have column names (species names) and row names (sampling plots)
##
## The function c.score calculates the C-score matrix to detect species association, for the whole community and for species pairs
## Randomization: column sum is fixed

## It returns the C-score index for the observed community (ObsCscoreTot), p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the working directory where the same 
## metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version)

## NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli and McCabe 2002 - Ecology)
##
## Literature 
## Gotelli, N.J. and D.J. McCabe. 2002. Ecology, 83, 2091-2096.
## Gotelli, N.J. and Ulrich, W. 2010. Oecologia, 162, 463-477
## Stone, L. & Roberts, A. 1990. Oecologia, 85, 74-79
# library(vegan); library(ade4) ;dependencies: permatswap::vegan; randtest::ade4


ecospat.Cscore <- function(data, nperm, outpath)
{

  # C-coef Observed matrix
  cat("Computing observed co-occurence matrix", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)


  spec.occ <- data.matrix(data)

  #### C-score calculation
  coocc <- t(spec.occ) %*% spec.occ
  n.spec <- dim(coocc)[1]
  mat1 <- array(apply(spec.occ, MARGIN = 2, sum), dim = c(n.spec, n.spec))
  mat2 <- t(array(apply(spec.occ, MARGIN = 2, sum), dim = c(n.spec, n.spec)))
  mat.obs.c.coef <- (mat1 - coocc) * (mat2 - coocc)  # observed c score
  df.obs.c.coef <- data.frame(Col = rep(1:ncol(mat.obs.c.coef), each = ncol(mat.obs.c.coef)),
    Row = rep(1:nrow(mat.obs.c.coef), nrow(mat.obs.c.coef)), Sp1 = rep(colnames(mat.obs.c.coef),
      each = ncol(mat.obs.c.coef)), Sp2 = rep(rownames(mat.obs.c.coef), nrow(mat.obs.c.coef)),
    Co.Occ = c(mat.obs.c.coef))
  v.diago.inf <- c(rownames(df.obs.c.coef)[df.obs.c.coef[, 1] > df.obs.c.coef[,
    2]], rownames(df.obs.c.coef)[df.obs.c.coef[, 1] == df.obs.c.coef[, 2]])
  df.obs.c.coef <- df.obs.c.coef[-as.numeric(v.diago.inf), ]
  CscoreTot <- mean(df.obs.c.coef$Co.Occ)

  ### Matrix to store the permutations
  mat.perm <- matrix(0, nrow(df.obs.c.coef), nperm, dimnames = list(c(paste(df.obs.c.coef[,
    3], df.obs.c.coef[, 4])), c(1:nperm)))


  ### Permutations C-score

  cat("Computing permutations", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)

  for (i in 1:nperm) {
    if (i == 1) {
      cat(nperm, " permutations to go", "\n", append = FALSE)
      cat(".............", "\n", append = FALSE)
    }
    if (i == nperm/2) {
      cat(nperm/2, " permutations to go", "\n", append = FALSE)
      cat(".............", "\n", append = FALSE)
    }


    spec.occ.perm1 <- data.matrix(data)
    spec.occ.perm1 <- permatswap(spec.occ.perm1, fixedmar = "both", mtype = "prab",
      times = 1)  # times=1 : separate swapping sequence that always begins with the original matrix
    spec.occ.perm <- as.matrix(spec.occ.perm1[[3]][[1]])

    coocc.perm <- t(spec.occ.perm) %*% spec.occ.perm
    mat1.perm <- array(apply(spec.occ.perm, MARGIN = 2, sum), dim = c(n.spec,
      n.spec))
    mat2.perm <- t(array(apply(spec.occ.perm, MARGIN = 2, sum), dim = c(n.spec,
      n.spec)))

    mat.obs.c.coef.perm <- (mat1.perm - coocc.perm) * (mat2.perm - coocc.perm)

    df.obs.c.coef.perm <- data.frame(Col = rep(1:ncol(mat.obs.c.coef.perm), each = ncol(mat.obs.c.coef.perm)),
      Row = rep(1:nrow(mat.obs.c.coef.perm), nrow(mat.obs.c.coef.perm)), Sp1 = rep(colnames(mat.obs.c.coef),
        each = ncol(mat.obs.c.coef.perm)), Sp2 = rep(rownames(mat.obs.c.coef),
        nrow(mat.obs.c.coef.perm)), Co.Occ = c(mat.obs.c.coef.perm))



    df.obs.c.coef.perm <- df.obs.c.coef.perm[-as.numeric(v.diago.inf), ]

    # Store result of permuation
    mat.perm[, i] <- df.obs.c.coef.perm[, 5]

  }


  #### Calculate C-score and SES for the whole community

  vec.CScore.tot <- as.vector(apply(mat.perm, MARGIN = 2, mean))
  SimulatedCscore <- mean(vec.CScore.tot)
  sd.SimulatedCscore <- sd(vec.CScore.tot)
  ses <- (CscoreTot - SimulatedCscore)/sd.SimulatedCscore
  randtest.less <- as.randtest(vec.CScore.tot, CscoreTot, alter = "less")
  pval.less <- randtest.less$pvalue
  randtest.greater <- as.randtest(vec.CScore.tot, CscoreTot, alter = "greater")
  pval.greater <- randtest.greater$pvalue
  # plot(randtest.greater, xlab= 'Simulated C-scores',main=paste('', sep=''))


  ### Calculate P-values based on random distribution

  mat.pval <- matrix(0, nrow(mat.perm), 4, dimnames = list(rownames(mat.perm), c("Obs.Co.Occ",
    "SES_Cscore", "pval_less", "pval_greater")))
  mat.pval[, 1] <- df.obs.c.coef[, 5]

  cat("Computing P-values", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)

  for (k in 1:nrow(mat.perm)) {

    mat.pval[k, 2] <- (df.obs.c.coef[k, 5] - mean(mat.perm[k, ]))/sd(mat.perm[k,
      ])
    randtest <- as.randtest(sim = mat.perm[k, ], obs = df.obs.c.coef[k, 5], alter = "less")
    mat.pval[k, 3] <- randtest$pvalue
    randtest <- as.randtest(sim = mat.perm[k, ], obs = df.obs.c.coef[k, 5], alter = "greater")
    mat.pval[k, 4] <- randtest$pvalue
  }


  ### Exporting Co-occ matrix
  cat("Exporting dataset", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)

  hist(as.vector(mat.pval[, 2]), xlab = "ses", main = paste("Histogram of standardized effect size"))
  abline(v = c(2, -2), col = "red")

  mat.pval.names <- data.frame(df.obs.c.coef[, 3:4], mat.pval, df.obs.c.coef.perm[,
    5])
  mat.pval.names2 <- data.frame(mat.pval.names[, 1:3], mat.pval.names[, 7], mat.pval.names[,
    4:6])
  names(mat.pval.names2)[3] <- "obs.C-score"
  names(mat.pval.names2)[4] <- "exp.C-score"
  write.table(mat.pval.names2, file = paste(outpath, "/Cscores.txt", sep = ""),
    sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE)

  tab <- mat.pval.names2
  v <- c(0)
  for (i in 1:nrow(tab)) {
    if (tab[i, 6] <= 0.05 || tab[i, 7] <= 0.05) {
      v <- c(v, i)
    }
  }
  m <- data.frame()
  for (j in 1:length(v)) {
    m <- rbind(m, tab[v[j], ])
  }

  m1 <- na.omit(m)

  write.table(m1, file = paste(outpath, "/Sign_Cscores.txt", sep = ""), sep = "\t",
    append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE)
  l <- list(ObsCscoreTot = CscoreTot, SimCscoreTot = SimulatedCscore, PVal.less = pval.less,
    PVal.greater = pval.greater, SES.Tot = ses)

  return(l)

  cat("End computation", "\n", append = FALSE)

  # END FUNCTION
}


#########        	CO-OCCURRENCE ANALYSIS			     ########
#########  					         &         					   ########
#########  ENVIRONMENTALLY CONSTRAINED NULL MODELS ########
#########                  CtRA1                   ########
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
# p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the path specified where the same 
# metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version)
# # NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002 - Ecology)
#
# presence: presence absence table
# pred: table with values of probability of presence for each species in each site
# nperm: number of permutation in the null model
# outpath: path to specify where to save the tables
# NB: Format required for imput databases: a plots (rows) x species (columns) matrix
# Input matrices should have column names (species names) and row names (sampling plots)


##################################################################################################################################################

# NullModels(presence, pred, 1000, outpath) # launch command
# library(ade4)

#############################
## Constrained permutation ##
#############################
Constrperm <- function(presence, pred) {
  nbsps <- ncol(presence)  # species number
  nbsites <- nrow(presence)  # sites number
  
  nbocc <- as.vector(apply(presence, MARGIN = 2, sum))  # occurrences number for each species
  sumprob <- as.vector(apply(pred, MARGIN = 2, sum))  # occurrence probability sum for each species
  
  matsum <- matrix(0, ncol = nbsps, nrow = nbsites)
  
  for (i in 1:nbsps) {
    matsum[, i] <- sumprob[i]
  }
  
  # creates a relative probability matrix
  pred <- pred/matsum
  
  transpo <- t(as.matrix(presence))
  sps.names <- row.names(transpo)
  sites <- row.names(presence)
  noms <- list(sites, sps.names)
  
  # creates a random matrix weighted by the species' probabilities of occurrence
  nullmat <- matrix(0, nrow = nbsites, ncol = nbsps, dimnames = noms)
  randr <- matrix(runif(nbsites * nbsps), nbsites)
  randr <- randr * pred
  
  for (i in 1:nbsps) {
    a <- as.vector(randr[, i])
    names(a) <- (1:nbsites)
    a <- sort(a)
    b <- as.numeric(names(a))
    x <- nbsites - nbocc[i] + 1
    
    for (j in x:nbsites) {
      r <- b[j]
      nullmat[r, i] <- 1
    }
  }
  return(nullmat)
}


SpeciesCooccurrenceStats <- function(presence) {
  nbsps <- ncol(presence)
  nbsites <- nrow(presence)
  nbocc <- as.vector(apply(presence, MARGIN = 2, sum))
  
  presence <- as.matrix(presence)
  coocc <- t(presence) %*% presence
  nbspec <- dim(coocc)[1]
  mat1 <- array(apply(presence, MARGIN = 2, sum), dim = c(nbsps, nbsps))
  mat2 <- t(array(apply(presence, MARGIN = 2, sum), dim = c(nbsps, nbsps)))
  Cscoreperspecies <- (mat1 - coocc) * (mat2 - coocc)
  
  # C-score on all matrix
  df.synthesis1 <- data.frame(Col = rep(1:ncol(Cscoreperspecies), each = ncol(Cscoreperspecies)), 
    Row = rep(1:nrow(Cscoreperspecies), nrow(Cscoreperspecies)), Sps1 = rep(colnames(Cscoreperspecies), 
      each = ncol(Cscoreperspecies)), Sps2 = rep(rownames(Cscoreperspecies), 
      nrow(Cscoreperspecies)), CScore = c(Cscoreperspecies))
  v.diago.inf <- c(rownames(df.synthesis1)[df.synthesis1[, 1] > df.synthesis1[, 
    2]], rownames(df.synthesis1)[df.synthesis1[, 1] == df.synthesis1[, 2]])
  df.synthesis <- df.synthesis1[-as.numeric(v.diago.inf), ]
  Cscore <- mean(df.synthesis[, 5])
  
  return(list(coocc = Cscoreperspecies, synth = df.synthesis1, synth_lt = df.synthesis, 
    Cscore = Cscore))
  
}  ### END function ###


################# Null models ##

ecospat.cons_Cscore <- function(presence, pred, nperm, outpath) {
  cat("Computing observed co-occurence matrix", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  
  nbsps <- ncol(presence)
  nbsites <- nrow(presence)
  
  # Co-occurrence statistics on observed matrix
  Obs <- SpeciesCooccurrenceStats(presence)
  
  CooccObs <- Obs$coocc
  synthObs <- Obs$synth_lt  #summary table with Co-occurrence values for each couple
  CscoreTot <- Obs$Cscore
  
  CooccProb <- matrix(0, nrow = nrow(synthObs), ncol = nperm, dimnames = list(c(paste(synthObs[, 
    3], synthObs[, 4])), c(1:nperm)))
  
  cat("Computing permutations", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  
  if (nperm > 0) {
    for (z in 1:nperm) {
      # cat(z, '\n',append =TRUE)
      degenerated <- TRUE
      while (degenerated == TRUE) {
        # looking for a matrix without degenerated sites
        random.distrib.matrix <- Constrperm(presence, pred)
        sumSites <- as.vector(apply(random.distrib.matrix, MARGIN = 2, sum))
        for (i in 1:nbsites) {
          ifelse(sumSites[i] == 0, degenerated <- TRUE, degenerated <- FALSE)
        }
      }
      
      # Co-occurrence statistics on randomized matrix
      Rnd <- SpeciesCooccurrenceStats(random.distrib.matrix)
      
      synthRnd <- Rnd$synth_lt
      
      CooccProb[, z] <- synthRnd[, 5]  # pairwise C-score for simulated matrix 
    }
    
    ## for the whole community
    vec.CScore.tot <- as.vector(apply(CooccProb, MARGIN = 2, mean))  # C-score for all null communities (mean on the columns)
    SimulatedCscore <- mean(vec.CScore.tot)  # mean of Simulation C-score: Simulated C-score
    sd.SimulatedCscore <- sd(vec.CScore.tot)  # standard deviation of null communities
    ses <- (CscoreTot - SimulatedCscore)/sd.SimulatedCscore  # standardized effect size: It scales the results in units of standard deviations, 
    # which allows for meaningful comparisons among different tests
    
    randtest.less <- as.randtest(vec.CScore.tot, CscoreTot, alter = "less")
    pval.less <- randtest.less$pvalue
    randtest.greater <- as.randtest(vec.CScore.tot, CscoreTot, alter = "greater")
    pval.greater <- randtest.greater$pvalue
    plot(randtest.greater, xlab = "Simulated C-scores", main = paste("", sep = ""))
    
    ## for species pairs
    vec.Cscore.pairs <- as.vector(apply(CooccProb, MARGIN = 1, mean))  # Mean of null communities C-score for any pair of species
    mat_pval <- matrix(0, nrow(synthRnd), 5)
    mat_pval[, 1] <- synthObs[, 5]  # Observed C-score in the first column
    mat_pval[, 2] <- vec.Cscore.pairs  # Mean of simulated C-scores in the second column
    
    for (i in 1:nrow(CooccProb)) {
      randtest.less <- as.randtest(CooccProb[i, ], mat_pval[i, 1], alter = "less")
      mat_pval[i, 3] <- randtest.less$pvalue
      randtest.greater <- as.randtest(CooccProb[i, ], mat_pval[i, 1], alter = "greater")
      mat_pval[i, 4] <- randtest.greater$pvalue
      mat_pval[i, 5] <- (mat_pval[i, 1] - mean(CooccProb[i, ]))/(sd(CooccProb[i, 
        ]))  #ses for any pair of species
    }
  }
  
  cat("Permutations finished", date(), "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat("Exporting dataset", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  cat(".............", "\n", append = FALSE)
  
  df.synthesis <- data.frame(Col = rep(1:ncol(CooccObs), each = ncol(CooccObs)), 
    Row = rep(1:nrow(CooccObs), nrow(CooccObs)), Sps1 = rep(colnames(CooccObs), 
      each = ncol(CooccObs)), Sps2 = rep(rownames(CooccObs), nrow(CooccObs)), 
    Cooc_score = c(CooccObs))
  v.diago.inf <- c(rownames(df.synthesis)[df.synthesis[, 1] > df.synthesis[, 2]], 
    rownames(df.synthesis)[df.synthesis[, 1] == df.synthesis[, 2]])
  df.synthesis <- df.synthesis[-as.numeric(v.diago.inf), ]
  df.synthesis <- cbind(df.synthesis, mat_pval[, 2:5])
  names(df.synthesis)[5:9] <- c("C-scoreObs", "C-scoreExp", "p.less", "p.greater", 
    "ses")
  
  hist(as.vector(df.synthesis[, 9]), xlab = "ses", main = paste("Histogram of standardized effect size"))
  abline(v = c(2, -2), col = "red")
  # write.table(df.synthesis,file=paste(outpath,'Results_const_Cscores.txt',
  # sep=''),sep='\t',append= FALSE,row.names= FALSE,col.names=TRUE,quote=FALSE)
  
  # selection of rows with <=0.05
  tab <- df.synthesis
  v <- c(0)
  for (i in 1:nrow(tab)) {
    if (tab[i, 7] <= 0.05 || tab[i, 8] <= 0.05) {
      v <- c(v, i)
    }
  }
  m <- data.frame()
  for (j in 1:length(v)) {
    m <- rbind(m, tab[v[j], ])
  }
  
  m1 <- na.omit(m)
  
  write.table(m1, file = paste(outpath, "/Signific_const_Cscores.txt", sep = ""), 
    sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  l <- list(ObsCscoreTot = CscoreTot, SimCscoreTot = SimulatedCscore, PVal.less = pval.less, 
    PVal.greater = pval.greater, SES.Tot = ses)
  return(l)
  
}