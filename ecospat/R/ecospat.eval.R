################ MEVA.TABLE: Function originally from A. Guisan (Unil-ECOSPAT, Switzerland)
################ Modified from L. Maiorano and O. Broennimann (Unil-ECOSPAT, Switzerland)
ecospat.meva.table <- function(Pred, Sp.occ, th) # Pred: vector of predicted probabilities Sp.occ: vector of binary observations
# th: threshold used to cut the probability to binary values
{
  # eva <- list()
  a <- table(Pred >= th, Sp.occ)[4]
  b <- table(Pred >= th, Sp.occ)[2]
  c <- table(Pred >= th, Sp.occ)[3]
  d <- table(Pred >= th, Sp.occ)[1]
  N <- a + b + c + d
  tab <- table(Pred >= th, Sp.occ, dnn = list("Predicted values", "Observed values"))
  prev <- (a + c)/N
  ccr <- (a + d)/N
  mcr <- (b + c)/N
  se <- a/(a + c)
  sp <- d/(b + d)
  fpr <- b/(b + d)
  fnr <- c/(a + c)
  ppp <- a/(a + b)
  npp <- d/(c + d)
  or <- (a * d)/(c * b)
  kappa <- ecospat.cohen.kappa(tab)$kap
  nmi <- (-(a * log(a)) - (b * log(b)) - (c * log(c)) - (d * log(d)) + ((a + b) *
    log(a + b)) + ((c + d) * log(c + d)))/((N * log(N)) - ((a + c) * log(a +
    c)) - ((b + d) * log(b + d)))
  tss <- se + sp - 1
  mtrcNAME <- c("Prevalence", "Correct classification rate", "Mis-classification rate",
    "Sensitivity", "Specificity", "Positive predictive power", "Negative predictive power",
    "False positive rate", "False negative rate", "Odds Ratio", "Kappa", "Normalized mutual information",
    "True skill statistic")
  mtrcVALUE <- c(round(prev, digits = 4), round(ccr, digits = 4), round(mcr, digits = 4),
    round(se, digits = 4), round(sp, digits = 4), round(fpr, digits = 4), round(fnr,
      digits = 4), round(ppp, digits = 4), round(npp, digits = 4), round(or,
      digits = 4), round(kappa, digits = 4), round(nmi, digits = 4), round(tss,
      digits = 4))
  mtrcMATR <- matrix(c(mtrcNAME, mtrcVALUE), nrow = length(mtrcVALUE), ncol = 2,
    dimnames = list(c(1:13), c("Metric", "Value")))
  return(list(CONTINGENCY_TABLE = tab, EVALUATION_METRICS = mtrcMATR))
}


############### MAX-KAPPA ## Function originally from A. Guisan (Unil-ECOSPAT, Switzerland)
############### Modified from L. Maiorano (Unil-ECOSPAT, Switzerland) and Olivier Broennimann
ecospat.max.kappa <- function(Pred, Sp.occ) # Pred: vector of predicted probabilities Sp.occ: vector of binary observations
{
  FUN<-function(thresh){
    xtab<-table(Pred >= thresh, Sp.occ)
    return(ecospat.cohen.kappa(xtab)$kap )
  }
  
  threshold<-(1:100)/100
  k<-sapply(threshold,FUN)
  table<-data.frame(cbind(threshold,k))
  
  max.Kappa <- max(table$k,na.rm = TRUE)
  max.threshold <- table$threshold[which.max(table$k)][1]
  return(list(table=table,max.Kappa=max.Kappa,max.threshold=max.threshold))
}


############### MAX-TSS ## Function originally from L. Maiorano (Unil-ECOSPAT, Switzerland)
############### modyfying MAX-KAPPA from A. Guisan (Unil-ECOSPAT, Switzerland), modified by Olivier Broennimann
ecospat.max.tss <- function(Pred, Sp.occ) # Pred: vector of predicted probabilities Sp.occ: vector of binary observations
{
  FUN<-function(thresh){
    xtab<-table(Pred >= thresh, Sp.occ)
    a <- xtab[4]
    b <- xtab[2]
    c <- xtab[3]
    d <- xtab[1]
    se <- a/(a + c)
    sp <- d/(b + d)
    return(se + sp - 1)
  }
  threshold<-(1:100)/100
  tss<-sapply(threshold,FUN)
  table<-data.frame(cbind(threshold,tss))
  
  max.TSS <- max(table$tss,na.rm = TRUE)
  max.threshold <- table$threshold[which.max(table$tss)][1]
  return(list(table=table,max.TSS=max.TSS,max.threshold=max.threshold))
}

################## PLOT-K ## Function originally from L. Maiorano (Unil-ECOSPAT, Switzerland)
ecospat.plot.kappa <- function(Pred, Sp.occ) # Pred: vector of predicted probabilities Sp.occ: vector of binary observations
{
  k <- 0.01
  i <- 1
  evak <- data.frame(k)
  while (k < 1) {
    a <- table(Pred >= k, Sp.occ)[4]
    b <- table(Pred >= k, Sp.occ)[2]
    c <- table(Pred >= k, Sp.occ)[3]
    d <- table(Pred >= k, Sp.occ)[1]
    N <- a + b + c + d
    tab <- table(Pred >= k, Sp.occ)
    evak[i, "THRESHOLD"] <- k
    evak[i, "k"] <- ecospat.cohen.kappa(tab)$kap
    # evak[i,'k'] <- ecospat.cohen.kappa(table(Pred >= k, Sp.occ))$kap
    k <- k + 0.01
    i <- i + 1
  }
  plot(evak$THRESHOLD, evak$k, type = "l", xlab = "THRESHOLD", ylab = "KAPPA")
}


################## PLOT-TSS ## Function originally from L. Maiorano (Unil-ECOSPAT, Switzerland)
ecospat.plot.tss <- function(Pred, Sp.occ) # Pred: vector of predicted probabilities Sp.occ: vector of binary observations
{
  tss <- 0.01
  i <- 1
  evatss <- data.frame(tss)
  while (tss < 1) {
    a <- table(Pred >= tss, Sp.occ)[4]
    b <- table(Pred >= tss, Sp.occ)[2]
    c <- table(Pred >= tss, Sp.occ)[3]
    d <- table(Pred >= tss, Sp.occ)[1]
    N <- a + b + c + d
    se <- a/(a + c)
    sp <- d/(b + d)
    evatss[i, "THRESHOLD"] <- tss
    evatss[i, "tss"] <- se + sp - 1
    tss <- tss + 0.01
    i <- i + 1
  }
  plot(evatss$THRESHOLD, evatss$tss, type = "l", xlab = "THRESHOLD", ylab = "TRUE SKILL STATISTIC")
}

############################################################################################################################


ecospat.cohen.kappa <- function(xtab) {
  # computes Cohen's kappa and variance estimates, 95% CI from a symmetric
  # agreement table xtab return value is a list with elements kap and vark See
  # Bishop et al. Disc. Mult. Analysis pp. 395-397 1975

  # Originally coded by C. Randin, UNIL Adjusted by N.E.Zimmermann, WSL

  ci <- vector(length = 2)
  kap <- 0
  vark <- 0
  totn <- 0
  if (nrow(xtab) != ncol(xtab)) {
    # cat('\n freq. table not symmetric.\n')
    result = list(kap = kap, vark = vark, totn = totn)
    return(result)
  }
  # compute obs props:
  totn <- sum(xtab)
  # calc marginals, expected diag props.
  obs <- xtab/totn
  pi <- apply(xtab, 1, sum)/totn
  pj <- apply(xtab, 2, sum)/totn
  exp <- outer(pi, pj)
  theta.2 <- sum(diag(exp))
  theta.1 <- sum(diag(obs))
  kap <- (theta.1 - theta.2)/(1 - theta.2)
  theta.3 <- diag(obs) %*% (pi + pj)
  theta.4 <- sum(obs * outer(pi, pj, FUN = "+")^2)
  # cat('\n thetas 1,2,3,4 \n', theta.1, theta.2, theta.3, theta.4, '\n') calc
  # the var est
  vark <- (theta.1 * (1 - theta.1))/(1 - theta.2)^2
  vark <- vark + (2 * (1 - theta.1) * (2 * theta.1 * theta.2 - theta.3))/(1 - theta.2)^3
  vark <- vark + ((1 - theta.1)^2 * (theta.4 - 4 * theta.2^2))/(1 - theta.2)^4
  vark <- as.vector(vark/totn)
  ci[1] <- kap - 1.96 * sqrt(vark)
  ci[2] <- kap + 1.96 * sqrt(vark)
  result = list(kap = kap, vark = vark, totn = totn, ci = ci)
  return(result)
}
