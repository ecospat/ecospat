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
library(dismo)
data('ecospat.testData')

# data for Soldanella alpina
data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")]

# copy maxent.jar file in the right folder
path.from<-system.file("extdata", "maxent.txt", package="ecospat")
path.to <- paste0(system.file(package="dismo"), "/java/maxent.txt")
path.to.renamed <- paste0(system.file(package="dismo"), "/java/maxent.jar")
file.copy(path.from,path.to,overwrite = TRUE)
file.rename(path.to, path.to.renamed)

if (file.exists(path.to.renamed) & require(rJava)) {
  me <- maxent(data.Solalp[,-1],data.Solalp[,1])
  ecospat.maxentvarimport (model=me, dfvar=data.Solalp[,-1], nperm=5)
  }
}
