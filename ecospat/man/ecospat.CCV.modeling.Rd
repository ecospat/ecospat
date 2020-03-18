\name{ecospat.CCV.modeling}
\alias{ecospat.CCV.modeling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs indivudual species distribuion models with SDMs or ESMs
%%  ~~function to do ... ~~
}
\description{
Creates probabilistic prediction for all species based on SDMs or ESMs and returns their evaluation metrics and variable importances.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ecospat.CCV.modeling(sp.data, 
                     env.data, 
                     xy,
                     DataSplitTable=NULL,
                     DataSplit = 70, 
                     NbRunEval = 25,
                     minNbPredictors =5,
                     validation.method = "cross-validation",
                     models.sdm = c("GLM","RF"), 
                     models.esm = "CTA", 
                     modeling.options.sdm = NULL, 
                     modeling.options.esm = NULL, 
                     ensemble.metric = "AUC", 
                     ESM = "YES",
                     parallel = FALSE, 
                     cpus = 4,
                     VarImport = 10,
                     modeling.id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.data}{a data.frame where the rows are sites and the columns are species (values 1,0)}
  \item{env.data}{either a data.frame where rows are sites and colums are environmental variables or a raster stack of the envrionmental variables}
  \item{xy}{two column data.frame with X and Y coordinates of the sites (most be same coordinate system as \code{env.data})}
  \item{DataSplitTable}{a table providing \code{TRUE/FALSE} to indicate what points are used for calibration and evaluation. As returned by \code{\link{ecospat.CCV.createDataSplitTable}}}
  \item{DataSplit}{percentage of dataset observations retained for the model training  (only needed if no \code{DataSplitTable} provided)}
  \item{NbRunEval}{number of cross-validatio/split sample runs (only needed if no \code{DataSplitTable} provided)}
  \item{minNbPredictors}{minimum number of occurences [min(presences/Absences] per predicotors needed to calibrate the models}
  \item{validation.method}{either "cross-validation" or "split-sample" used to validate the communtiy predictions (only needed if no \code{DataSplitTable} provided)}
  \item{models.sdm}{modeling techniques used for the normal SDMs. Vector of models names choosen among \code{'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips' and 'MAXENT.Tsuruoka'}}
  \item{models.esm}{modeling techniques used for the ESMs. Vector of models names choosen among \code{'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips' and 'MAXENT.Tsuruoka'}}
  \item{modeling.options.sdm}{modeling options for the normal SDMs. \code{"BIOMOD.models.options"}" object returned by \code{\link[biomod2]{BIOMOD_ModelingOptions}}}
  \item{modeling.options.esm}{modeling options for the ESMs. \code{"BIOMOD.models.options"} object returned by \code{\link[biomod2]{BIOMOD_ModelingOptions}}}
  \item{ensemble.metric}{evaluation score used to weight single models to build ensembles: \code{'AUC', 'Kappa' or 'TSS'}}
  \item{ESM}{either \code{'YES'} (ESMs allowed), \code{'NO'} (ESMs not allowed) or \code{'ALL'} (ESMs used in any case)}
  \item{parallel}{should parallel computing be allowed (\code{TRUE/FALSE})}
  \item{cpus}{number of cpus to use in parallel computing}
  \item{VarImport}{number of permutation runs to evaluate variable importance}
  \item{modeling.id}{character, the ID (=name) of modeling procedure. A random number by default}
}
\details{
The basic idea of the community cross-validation (CCV) is to use the same data (sites) for the model calibration/evaluation of all species. This ensures that there is "independent" cross-validation/split-sample data available not only at the individual species level but also at the community level. This is key to allow an unbiased estimation of the ability to predict species assemblages (Scherrer et al. 2018). 
The output of the ecospat.CCV.modeling function can then be used to evaluate the species assemblage predictions with the \code{\link{ecospat.CCV.communityEvaluation.bin}} or \code{\link{ecospat.CCV.communityEvaluation.prob}} functions.
}
\value{
  \item{modelling.id}{character, the ID (=name) of modeling procedure}
  \item{output.files}{vector with the names of the files written to the hard drive}
  \item{speciesData.calibration}{a 3-dimensional array of presence/absence data of all species for the calibration plots used for each run}
  \item{speciesData.evaluation}{a 3-dimensional array of presence/absence data of all species for the evaluation plots used for each run}
  \item{speciesData.full}{a data.frame of presence/absence data of all species (same as \code{sp.data} input)}
  \item{DataSplitTable}{a matrix with \code{TRUE/FALSE} for each model run (\code{TRUE}=Calibration point, \code{FALSE}=Evaluation point)}
  \item{singleSpecies.ensembleEvaluationScore}{a 3-dimensional array of single species evaluation metrics (\code{'Max.KAPPA', 'Max.TSS', 'AUC of ROC'})}
  \item{singleSpecies.ensembleVariableImportance}{a 3-dimensional array of single species variable importance for all predictors}
  \item{singleSpecies.calibrationSites.ensemblePredictions}{a 3-dimensional array of the predictions for each species and run at the calibration sites}
  \item{singleSpecies.evaluationSites.ensemblePredictions}{a 3-dimensional array of the predictions for each species and run at the evaluation sites}
  \item{allSites.averagePredictions.cali}{a matrix with the average predicted probabilities for each site across all the runs the sites were used for model calibration}
  \item{allSites.averagePredictions.eval}{a matrix with the average predicted probabilities for each site across all the runs the sites were used as independent evaluation sites}
}
\references{
Scherrer, D., D'Amen, M., Mateo, M.R.G., Fernandes, R.F. & Guisan , A. (2018) How to best threshold and validate stacked species assemblages? Community optimisation might hold the answer. Methods in Ecology and Evolution, in review
}
\author{
Daniel Scherrer <daniel.j.a.scherrer@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ecospat.CCV.createDataSplitTable}}; \code{\link{ecospat.CCV.communityEvaluation.bin}}; \code{\link{ecospat.CCV.communityEvaluation.prob}}
}
\examples{
\donttest{
#Loading species occurence data and remove empty communities
testData <- ecospat.testData[,c(24,34,43,45,48,53,55:58,60:63,65:66,68:71)]
sp.data <- testData[which(rowSums(testData)>0), sort(colnames(testData))]

#Loading environmental data
env.data <- ecospat.testData[which(rowSums(testData)>0),4:8]

#Coordinates for all sites
xy <- ecospat.testData[which(rowSums(testData)>0),2:3]

#Running all the models for all species
myCCV.Models <- ecospat.CCV.modeling(sp.data = sp.data,
                                     env.data = env.data,
                                     xy = xy,
                                     NbRunEval = 5,
                                     minNbPredictors = 10,
                                     VarImport = 3)
}
}
