\name{ecospat.SESAM.prr}
\alias{ecospat.SESAM.prr}

\title{SESAM Probability Ranking Rule
}
\description{
Implement the SESAM framework to predict community composition using a `probability ranking` rule.
}
\usage{ecospat.SESAM.prr(proba, sr, verbose = FALSE)
}
\arguments{
    \item{proba}{A data frame object of SDMs probabilities (or other sources) for all species. Column names (species names SDM) and row name (sampling sites) (need to have defined row names).}
    \item{sr}{
A data frame object with species richness value in the first column. Sites should be arranged in the same order as in the `prob` argument. 
    }
    \item{verbose}{Boolean indicating whether to print progress output during calculation. Default is FALSE.}
}
\details{
The SESAM framework implemented in ecospat is based on 1) probabilities of individual species presence for each site - these can be obtained for example by fitting SDMs. This step represents the application of an environmental filter to the community assembly, 2) richness predictions for each site - the richness prediction can be derived in different ways, for instance by summing probabilities from individual species presence for each site  or by fitting direct richness models. This step represents the application of a macroecological constraint to the number of species that can coexist in the considered unit, 3) a biotic rule to decide which species potentially present in the site are retained in the final prediction to match the richness value predicted. The biotic rule applied here is called `probability ranking` rule: the community composition in each site is determined by ranking the species in decreasing order of their predicted probability of presence from SDMs up to a richness prediction.
}
\value{
Returns a `.txt` file saved in the working directory that contains the community prediction by the SESAM framework, i.e. binary predictions for all species (columns) for each site (rows).
}
\author{
Manuela D`Amen \email{manuela.damen@unil.ch} and Anne Dubuis \email{anne.dubuis@gmail.com}
}
\references{
D`Amen, M., A. Dubuis, R.F. Fernandes, J. Pottier, L. Pellissier and A. Guisan. 2015. Using species richness and functional traits predictions to constrain assemblage predictions from stacked species distribution models. \emph{J. Biogeogr.}, \bold{42}, 1255-1266.

Guisan, A. and C. Rahbek. 2011. SESAM - a new framework integrating macroecological and species distribution models for predicting spatio-temporal patterns of species assemblages. \emph{J. Biogeogr.}, \bold{38}, 1433-1444.
}


\examples{
proba <- ecospat.testData[,73:92]
sr <- as.data.frame(rowSums(proba))
ppr<-ecospat.SESAM.prr(proba, sr)
head(ppr)

}
