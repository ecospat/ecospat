\name{ecospat.ESM.VarContrib}
\alias{ecospat.ESM.VarContrib}

\title{
Variable contribution in ESM
}
\description{
calculates the variable contribution of each variable and method in an ESM model
}
\usage{
    ecospat.ESM.VarContrib(ESM.modeling.output,
                           ESM_EF.output)
}

\arguments{
    \item{ESM.modeling.output}{\code{BIOMOD.formated.data} object returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
    \item{ESM_EF.output}{\code{BIOMOD.formated.data} object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
}
\details{
Calculates the ration between sum of weights of bivariate models where a focal variable was used and sum of weights of all bivariate models. This gives an indication on the proportional contribution of the variable in the final ensemble model.
In the case of multiple methods (e.g., GLM, GAM...), the contributions are counted per method. For ensemble model, the contributions are then weighted means (based on the weighting score as chosen in ecospat.ESM.EnsembleModeling()) of single methods
    
}

\value{
Returns a dataframe with contribution values (i.e. proportional contribution) by variable and model

}
\author{
Olivier Broennimann <Olivier.Broennimann@unil.ch> with contributions of Heidi Mod \email{Heidi.mod1@gmail.com} and Daniel Scherrer \email{aniel.j.a.scherrer@gmail.com}
}

\seealso{
\code{\link[ecospat]{ecospat.ESM.Modeling}}, \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}, \code{\link[ecospat]{ecospat.ESM.Projection}}, \code{\link[ecospat]{ecospat.ESM.EnsembleProjection}}

}
\examples{
   \donttest{       
require(biomod2)
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
                                models=c('GLM','RF'),
                                NbRunEval=2,
                                DataSplit=70,
                                Prevalence=0.5,
                                weighting.score=c("AUC"),
                                parallel=FALSE)  

### Evaluation and average of simple bivariate models to ESMs
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)

## get the model performance of ESMs 
my.ESM_EF$ESM.evaluations

## get the weights of the single bivariate models used to build the ESMs
my.ESM_EF$weights

## get the variable contributions of ESMs
ecospat.ESM.VarContrib(my.ESM,my.ESM_EF)
}
}
