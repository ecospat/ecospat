\name{ecospat.ESM.threshold}
\alias{ecospat.ESM.threshold}

\title{
Thresholds for Ensamble of Small Models
}
\description{
This function evaluates the full model which is used for projections and provides thresholds to produce binary maps.

}
\usage{
    ecospat.ESM.threshold( ESM.EnsembleModeling.output, 
                           PEplot = FALSE)
}
\arguments{
    \item{ESM.EnsembleModeling.output}{a list object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
    \item{PEplot}{logical. Should the predicted to expected ratio along the suitability class from the boyce index be plotted. Default FALSE (see \code{\link[ecospat]{ecospat.boyce}}} 

}
\details{
This function provides evaluation scores of the full model (no split sampling) and thresholds which can be used to convert suitability maps into binary maps. Various thresholds are provided: TSS (where sensitivity and specificity are maximised), MPA 1.0 (where all presences are prdicted positive), MPA 0.95 (where 95\% of all presences are predicted positive), MPA 0.90 (where 90\% of all presences are predicted positive), Boyce.th.min (the lowest suitability value where the predicted/expected ratio is >1) and Boyce.th.max (the highest suitability value where the predicted/expected ratio is =1). 
}

\value{
A data.frame with evluation scores and thresholds.
}

\author{
Frank Breiner \email{frank.breiner@unil.ch}

}
\references{
Hirzel, Alexandre H., et al. Evaluating the ability of habitat suitability models to predict species presences. \emph{Ecological modelling}, \bold{199.2} (2006): 142-152.

Engler, Robin, Antoine Guisan, and Luca Rechsteiner. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. \emph{Journal of applied ecology}, \bold{41.2} (2004): 263-274.

Fielding, Alan H., and John F. Bell. A review of methods for the assessment of prediction errors in conservation presence/absence models." \emph{Environmental conservation}, \bold{24.1} (1997): 38-49.
          }

\seealso{
\code{\link[ecospat]{ecospat.ESM.EnsembleProjection}}, \code{\link[ecospat]{ecospat.boyce}}, \code{\link[ecospat]{ecospat.mpa}},\code{\link[ecospat]{ecospat.binary.model}}
}
\examples{
   \donttest{
library(biomod2)

# Loading test data
data(ecospat.testNiche.inv)
inv <- ecospat.testNiche.inv

# species occurrences
xy <- inv[,1:2]
sp_occ <- inv[11]

# env
current <- inv[3:10]

### Formating the data with the BIOMOD_FormatingData() function from the package biomod2

sp <- 1
myBiomodData <- BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
                                      expl.var = current,
                                      resp.xy = xy,
                                      resp.name = colnames(sp_occ)[sp])

### Calibration of simple bivariate models
my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
                                models=c('GLM','CTA'),
                                NbRunEval=2,
                                DataSplit=70,
                                weighting.score=c('AUC'),
                                parallel=FALSE)  

### Evaluation and average of simple bivariate models to ESMs
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c('SomersD'),threshold=0)


output.TH <- ecospat.ESM.threshold(my.ESM_EF,PEplot = TRUE)

}}

