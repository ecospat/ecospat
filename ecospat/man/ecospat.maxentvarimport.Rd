\name{ecospat.maxentvarimport}

\alias{ecospat.maxentvarimport}

\title{Maxent Variable Importance}

\description{Calculate the importance of variables for Maxent in the same way Biomod does, by randomly permuting each predictor variable independently, and computing the associated reduction in predictive performance.}

\usage{ecospat.maxentvarimport (model, dfvar, nperm)}

\arguments{
  \item{model}{The name of the maxent model.}
  \item{dfvar}{A dataframe object with the environmental variables.}
  \item{nperm}{The number of permutations in the randomization process. The default is 5.}
}

\details{The calculation is made as biomod2 "variables_importance" function. It's more or less base on the same principle than randomForest variables importance algorithm. The principle is to shuffle a single variable of the given data. Make model prediction with this 'shuffled' data.set. Then we compute a simple correlation (Pearson's by default) between references predictions and the 'shuffled' one. The return score is 1-cor(pred_ref,pred_shuffled). The highest the value, the more influence the variable has on the model. A value of this 0 assumes no influence of that variable on the model. Note that this technique does not account for interactions between the variables.}

\value{
 a \code{list} which contains a \code{data.frame} containing variables importance scores for each permutation run. }

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}


\examples{
\dontrun{
ecospat.cv.example()
model <- get ("me.Achillea_millefolium", envir=ecospat.env)
dfvar <- ecospat.testData[4:8]
nperm <- 5
ecospat.maxentvarimport (model, dfvar, nperm)
}
}
