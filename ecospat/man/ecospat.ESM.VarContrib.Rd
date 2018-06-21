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
                           ESM_EF.output,
                           scaling)
}

\arguments{
    \item{ESM.modeling.output}{\code{BIOMOD.formated.data} object returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
    \item{ESM_EF.output}{\code{BIOMOD.formated.data} object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
    \item{scaling}{rescaling method: 'plusminus1','plusminus1bymodel','01','01bymodel' or 'none'. See details. default is "plusminus1"}
}
\details{
Calculates the difference in bivariate model weights were a focal variable was used compared to all bivariate model weights. 
It gives an indication on the contribution of the variable in the final ensemble model. Without rescaling the contributions 
are in the scale of the bivariate model weights. They can be rescale between 0 and 1 (scaling='01') or multiplied to a maximum 
of 1 (or -1 if the biggest contribution is negative; scaling ='plusminus1'). Scaling methods with 'bymodel' do the same 
but within each modeling techiques.
}

\value{
A dataframe with contribution values by variable and model

}
\author{
Olivier Broennimann <Olivier.Broennimann@unil.ch> with contributions of Heidi Mod \email{Heidi.mod1@gmail.com} and Daniel Scherrer \email{aniel.j.a.scherrer@gmail.com}
}

\seealso{
\code{\link[ecospat]{ecospat.ESM.Modeling}}, \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}, \code{\link[ecospat]{ecospat.ESM.Projection}}, \code{\link[ecospat]{ecospat.ESM.EnsembleProjection}}

}
\examples{
   \dontrun{
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
ecospat.ESM.VarContrib(my.ESM,my.ESM_EF,scaling="01")
ecospat.ESM.VarContrib(my.ESM,my.ESM_EF,scaling="01bymodel")

}}
