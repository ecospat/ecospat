#### exdet_function
#### Coded by Blaise Petitpierre (bpetitpierre@gmail.com), 9th December 2015 after Mesgaran et al 2014 
#### it allows to assess climate analogy
#### between a projection extent (p) and a reference extent (ref, usend in general
#### as the background to calibrate SDMS) ref : a dataframe with the value of the
#### variables (i.e columns) for each point of the reference exent p : a dataframe
#### with the value of the variables (i.e columns) for each point of the projection
#### exent return a vector. Values below 0 are novel conditions at the univariate
#### level (similar to the MESS), values between 0 and 1 are analog and values above
#### 1 are novel covariate condtions. For more information, see the reference
#### References: Mesgaran, M.B., Cousens, R.D. & Webber, B.L. (2014) Here be
#### dragons: a tool for quantifying novelty due to covariate range and correlation
#### change when projecting species distribution models. Diversity & Distributions,
#### 20: 1147-1159, DOI: 10.1111/ddi.12209

ecospat.exdet <- function(ref, p) {
  p <- as.matrix(p)
  a <- apply(ref, 2, min)
  b <- apply(ref, 2, max)
  minref <- matrix(a, nrow = nrow(p), ncol = ncol(p), byrow = TRUE)
  maxref <- matrix(b, nrow = nrow(p), ncol = ncol(p), byrow = TRUE)
  
  nt1 <- rowSums(apply(array(data = c(p - minref, maxref - p, rep(0, nrow(p) * 
    ncol(p))), dim = c(dim(p), 3)), c(1, 2), min)/(maxref - minref))
  
  tokeep <- which(nt1 == 0)
  p <- matrix(p[tokeep, ], ncol = ncol(p))
  a <- apply(ref, 2, mean)
  b <- var(ref)
  
  mah.ref <- mahalanobis(x = ref, center = a, cov = b)
  mah.pro <- mahalanobis(x = p, center = a, cov = b)
  mah.max <- max(mah.ref[is.finite(mah.ref)])
  nt2 <- mah.pro/mah.max
  nt1[tokeep] <- nt2
  return(nt1)
}