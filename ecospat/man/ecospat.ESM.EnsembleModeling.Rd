\name{ecospat.ESM.EnsembleModeling}
\alias{ecospat.ESM.EnsembleModeling}

\title{
Ensamble of Small Models: Evaluates and Averages Simple Bivariate Models To ESMs
}
\description{
This function evaluates and averages simple bivariate models by weighted means to Ensemble Small Models as in Lomba et al. 2010 and Breiner et al. 2015.

}
\usage{
    ecospat.ESM.EnsembleModeling( ESM.modeling.output, 
                                  weighting.score, 
                                  threshold=NULL, 
                                  models)
}
\arguments{
    \item{ESM.modeling.output}{a \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
    
    \item{weighting.score}{an evaluation score used to weight single models to build ensembles:"AUC","TSS","Boyce","Kappa","SomersD"

#the evaluation methods used to evaluate ensemble models ( see \code{"\link[=BIOMOD_Modeling]{BIOMOD_Modeling}"}  \bold{models.eval.meth} section for more detailed informations )
}
    \item{threshold}{
threshold value of an evaluation score to select the bivariate model(s) included for building the ESMs
    }
    \item{models}{
vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF','MAXENT.Phillips', "MAXENT.Tsuruoka" (same as in \code{biomod2})

#a character vector (either 'all' or a sub-selection of model names) that defines the models kept for building the ensemble models (might be useful for removing some non-preferred models)
    }

}
\details{
The basic idea of ensemble of small models (ESMs) is to model a species distribution based on small, simple models, for example all possible bivariate models (i.e. models that contain only two predictors at a time out of a larger set of predictors), and then combine all possible bivariate models into an ensemble (Lomba et al. 2010; Breiner et al. 2015).

The ESM set of functions could be used to build ESMs using simple bivariate models which are averaged using weights based on model performances (e.g. AUC) according to Breiner et al. (2015). They provide full functionality of the approach described in Breiner et al. (2015).


}
\value{

species:          species name
ESM.fit:          data.frame of the predicted values for the data used to build the models.
ESM.evaluations:  data.frame with evaluations scores for the ESMs
ESM.predictions:  Returns the projections of ESMs for the selected single models and their ensemble


A \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"}. This object will be later given to \code{\link[ecospat]{ecospat.ESM.EnsembleProjection}} if you want to make some projections of this ensemble-models.

}
\author{
Frank Breiner \email{frank.breiner@wsl.ch} with contributions of 
Olivier Broennimann \email{olivier.broennimann@unil.ch}

}
\references{

Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.
Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.
}

\seealso{
\code{\link[ecospat]{ecospat.ESM.Modeling}}, \code{\link[ecospat]{ecospat.ESM.Projection}}, \code{\link[ecospat]{ecospat.ESM.EnsembleProjection}}

\code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link{BIOMOD_Projection}}

}
\examples{
   \dontrun{
# Loading test data for the niche dynamics analysis in the invaded range
inv <- ecospat.testNiche.inv

# species occurrences
xy <- inv[,1:2]
sp_occ <- inv[11]

# env
current <- inv[3:10]



### Formating the data with the BIOMOD_FormatingData() function form the package biomod2
setwd(path.wd)
t1 <- Sys.time()
sp <- 1
myBiomodData <- BIOMOD_FormatingData( resp.var = as.numeric(sp_occ[,sp]),
                                      expl.var = current,
                                      resp.xy = xy,
                                      resp.name = colnames(sp_occ)[sp])

myBiomodOption <- Print_Default_ModelingOptions()


### Calibration of simple bivariate models
my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
                                models=c('GLM','RF'),
                                models.options=myBiomodOption,
                                NbRunEval=2,
                                DataSplit=70,
                                weighting.score=c("AUC"),
                                parallel=FALSE)  


### Evaluation and average of simple bivariate models to ESMs
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,weighting.score=c("SomersD"),threshold=0)

### Projection of simple bivariate models into new space 
my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
                                            new.env=current)

### Projection of calibrated ESMs into new space 
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
                                                        ESM.EnsembleModeling.output=my.ESM_EF)

}}

