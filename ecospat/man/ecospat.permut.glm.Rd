\name{ecospat.permut.glm}

\alias{ecospat.permut.glm}

\title{GLM Permutation Function}

\description{A permutation function to get p-values on GLM coefficients and deviance.}

\usage{ecospat.permut.glm (glm.obj, nperm, verbose = FALSE)}

\arguments{
  \item{glm.obj}{Any calibrated GLM or GAM object with a binomial error distribution.}
  \item{nperm}{The number of permutations in the randomization process.}
  \item{verbose}{Boolean indicating whether to print progress output during calculation. Default is FALSE.}
}

\details{Rows of the response variable are permuted and a new GLM is calibrated as well as deviance, adjusted deviance and coefficients are calculated. These random parameters are compared to the true parameters in order to derive p-values.}

\value{Return p-values that are how the true parameters of the original model deviate from the disribution of the random parameters. A p-value of zero means that the true parameter is completely outside the random distribution.}

\references{Hastie, T., R. Tibshirani and J. Friedman. 2001. \emph{Elements of Statistical Learning; Data Mining, Inference, and Prediction}, Springer-Verlag, New York.

Legendre, P. and L. Legendre. 1998. \emph{Numerical Ecology}, 2nd English edition. Elsevier Science BV, Amsterdam.}

\author{Christophe Randin \email{christophe.randin@unibas.ch}, Antoine Guisan \email{antoine.guisan@unil.ch} and Trevor Hastie}

\examples{

\donttest{
library(rms)
data('ecospat.testData')

# data for Soldanella alpina
data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")] 

# gbm model for Soldanella alpina

glm.Solalp <- glm("Soldanella_alpina ~ pol(ddeg,2) + pol(mind,2) + pol(srad,2) + pol(slp,2) 
      + pol(topo,2)", data = data.Solalp, family = binomial)
                  
# p-values
ecospat.permut.glm (glm.Solalp, 1000)
}

}