ecospat.binary.model <- function (Pred, Threshold)
{
  Threshold.Com <- c(0, Threshold, 0, Threshold, maxValue(Pred), 1)
  Threshold.Com.b <- matrix(Threshold.Com, ncol = 3, byrow = T)
  Pred.binary <- reclassify(Pred, Threshold.Com.b)
  names(Pred.binary) <- names(Pred)
  return(Pred.binary)
}
