\name{ecospat.ESM.Modeling}
\alias{ecospat.ESM.Modeling}

\title{
Ensamble of Small Models: Calibration of Simple Bivariate Models
}
\description{
This function calibrates simple bivariate models as in Lomba et al. 2010 and Breiner et al. 2015.
}
\usage{
    ecospat.ESM.Modeling( data, 
                          NbRunEval, 
                          DataSplit, 
                          DataSplitTable, 
                          Prevalence,
                          weighting.score, 
                          models, 
                          tune,
                          modeling.id, 
                          models.options, 
                          which.biva, 
                          parallel, 
                          cleanup,
                          Yweights)
}

\arguments{
    \item{data}{\code{BIOMOD.formated.data} object returned by \code{\link[biomod2]{BIOMOD_FormatingData}}
    }
    \item{NbRunEval}{
number of dataset splitting replicates for the model evaluation (same as in \code{biomod2})
    }
    \item{DataSplit}{
percentage of dataset observations retained for the model training (same as in \code{biomod2})
    }
    \item{DataSplitTable}{
a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column corresponds to a 'RUN'. If filled, arguments NbRunEval and DataSplit  will be ignored.
    }
    \item{Prevalence}{
either NULL or a 0-1 numeric used to build 'weighted response weights'. In contrast to Biomod the default is 0.5 (weighting presences equally to the absences). If NULL each observation (presence or absence) has the same weight (independent of the number of presences and absences).
    }
    \item{weighting.score}{
evaluation score used to weight single models to build ensembles: 'AUC', 'SomersD' (2xAUC-1), 'Kappa', 'TSS' or 'Boyce'
    }
    \item{models}{
vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF','MAXENT.Phillips', 'MAXENT.Tsuruoka' (same as in \code{biomod2})
    }
    \item{tune}{
logical. if true model tuning will be used to estimate optimal parameters for the models (Default: False).
    }    
    \item{modeling.id}{
character, the ID (=name) of modeling procedure. A random number by default.
    }
    \item{models.options}{
BIOMOD.models.options object returned by BIOMOD_ModelingOptions (same as in \code{biomod2})
    }
    \item{Yweights}{
response points weights. This argument will only affect models that allow case weights. 
    }
    \item{which.biva}{
integer. which bivariate combinations should be used for modeling? Default: all
    }
    \item{parallel}{
logical. If TRUE, the parallel computing is enabled (highly recommended)
    }
    \item{cleanup}{
numeric. Calls removeTmpFiles() to delete all files from rasterOptions()$tmpdir which are older than the given time (in hours). This might be necessary to prevent running over quota. No cleanup is used by default. 
    }
}
\details{
The basic idea of ensemble of small models (ESMs) is to model a species distribution based on small, simple models, for example all possible bivariate models (i.e. models that contain only two predictors at a time out of a larger set of predictors), and then combine all possible bivariate models into an ensemble (Lomba et al. 2010; Breiner et al. 2015).

The ESM set of functions could be used to build ESMs using simple bivariate models which are averaged using weights based on model performances (e.g. AUC) according to Breiner et al. (2015). They provide full functionality of the approach described in Breiner et al. (2015).

The argument \code{which.biva} allows to split model runs, e.g. if \code{which.biva} is 1:3, only the three first bivariate variable combinations will be modeled. This allows to run different biva splits on different computers. However, it is better not to use this option if all models are run on a single computer.
Default: running all biva models. 
NOTE: Make sure to give each of your biva runs a unique \code{modeling.id}. Please avoid space characters in your working directory path if you are using MAXENT.Phillips because this can cause an error.
}

\value{
A BIOMOD.models.out object (same as in \code{biomod2})
See \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} for details.

}
\author{
Frank Breiner \email{frank.breiner@wsl.ch} and Mirko Di Febbraro \email{mirkodifebbraro@gmail.com} with contributions of Olivier Broennimann \email{olivier.broennimann@unil.ch}
}
\references{
Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. 
\emph{Biological Conservation}, \bold{143},2647-2657.

Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.

Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
          }

\seealso{
\code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}, \code{\link[ecospat]{ecospat.ESM.Projection}}, \code{\link[ecospat]{ecospat.ESM.EnsembleProjection}}

\code{\link[biomod2]{BIOMOD_FormatingData}}, \code{\link[biomod2]{BIOMOD_ModelingOptions}}, \code{\link[biomod2]{BIOMOD_Modeling}},\code{\link[biomod2]{BIOMOD_Projection}}

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

## There are multiple packages that can be used for parallel computation (if parallel = TRUE). 
# One option is:
# library(doParallel);cl<-makeCluster(4);doParallel::registerDoParallel(cl) 

### Calibration of simple bivariate models
my.ESM <- ecospat.ESM.Modeling( data=myBiomodData,
                                models=c('GLM','RF'),
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

## get the model performance of ESMs 
my.ESM_EF$ESM.evaluations
## get the weights of the single bivariate models used to build the ESMs
my.ESM_EF$weights

## get the variable contributions of ESMs
ecospat.ESM.VarContrib(my.ESM,my.ESM_EF)
}}
