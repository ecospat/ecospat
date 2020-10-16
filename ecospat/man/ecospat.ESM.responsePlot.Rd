\name{ecospat.ESM.responsePlot}
\alias{ecospat.ESM.responsePlot}

\title{
Produce response plots for ESMs
}
\description{
This function creates response plots (evaluation strips) for Ensmebles of Small Models (ESMs).

}
\usage{
    ecospat.ESM.responsePlot( ESM.EnsembleModeling.output,
                              ESM.modeling.output,
                              fixed.var.metric = 'median')
}
\arguments{
    \item{ESM.modeling.output}{a list object returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
    \item{ESM.EnsembleModeling.output}{a list object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
    \item{fixed.var.metric}{either 'median' (default), 'mean', 'min' or 'max' specifying the statistic used to fix as constant the remaining variables when the predicted response is estimated for one of the variables. (same as in \code{\link[biomod2]{response.plot2}}} 
}

\details{
This function plots the response curves of a model for each variable, while keeping the remianing variables constant. This is an adaptation of the Evaluation Strip method proposed by Elith et al.(2005)
}

\value{
A plot of the response curves is produced (red line Ensemble, grey lines single algorithms) and a list with the output is provided.
}

\author{
Frank Breiner \email{frank.breiner@unil.ch}

}
\references{
Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting predicted responses from species distribution models. Ecological Modelling 186, 280-289.
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
output.plot <- ecospat.ESM.responsePlot(ESM.EnsembleModeling.output = my.ESM_EF,
                         ESM.modeling.output = my.ESM)

}}

