\name{ecospat.CCV.communityEvaluation.bin}
\alias{ecospat.CCV.communityEvaluation.bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates a range of community evaluation metrics based on different thresholding techniques.
}
\description{
The function uses the output of \code{\link{ecospat.CCV.modeling}} to calculate a range of community evaluation metrics based on a selection of thresholding techniques both for the calibration data and independent evaluation data.
}
\usage{
ecospat.CCV.communityEvaluation.bin(ccv.modeling.data,
                                    thresholds= c('MAX.KAPPA', 'MAX.ROC','PS_SDM'),
                                    community.metrics=c('SR.deviation','Sorensen'),
                                    parallel=FALSE,
                                    cpus=4,
                                    fix.threshold=0.5,
                                    MCE=5,
                                    MEM=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ccv.modeling.data}{a \code{'ccv.modeling.data'} object returned by \code{\link{ecospat.CCV.modeling}}}
  \item{thresholds}{a selection of thresholds (\code{'FIXED', 'MAX.KAPPA', 'MAX.ACCURACY', 'MAX.TSS', 'SENS_SPEC', 'MAX.ROC', 'OBS.PREVALENCE', 'AVG.PROBABILITY', 'MCE', 'PS_SDM, MEM'}) to be calculated and applied for the model evaluation.}
  \item{community.metrics}{a selection of community evaluation metrics (\code{'SR.deviation', 'community.AUC', 'community.overprediction' ,'community.underprediction',' community.accuracy', 'community.sensitivity', 'community.specificity', 'community.kappa', 'community.tss', 'Sorensen', 'Jaccard', 'Simpson'}) to be calculated for each seleted thresholding technique.}
  \item{parallel}{should parallel computing be allowed (\code{TRUE/FALSE})}
  \item{cpus}{number of cpus to use in parallel computing}
  \item{fix.threshold}{fixed threshold to be used. Only gets used if thresholding technique \code{FIXED} is selected.}
  \item{MCE}{maximum omission error (\%) allowed for the thresholding. Only gets used if thresholding technique \code{MCE} is selected.}
  \item{MEM}{a vetor with the species richness prediction of a MEM for each site. Only needed if \code{MEM} is selected.}
  
}
\details{
The function uses the probability output of the \code{\link{ecospat.CCV.modeling}} function and creates binary maps based on the selected thresholding methods. These binary maps are then used to calculate the selected community evaluation metrics both for the calibration and evaluation data of each modeling run.
}
\value{
  \item{DataSplitTable}{a matrix with \code{TRUE/FALSE} for each model run (\code{TRUE}=Calibration point, \code{FALSE}=Evaluation point)}
  \item{CommunityEvaluationMetrics.CalibrationSites}{a 4-dimensional array containing the community evaluation metrics for the calibartion sites of each run (\code{NA} means that the site was used for evaluation)}
  \item{CommunityEvaluationMetrics.EvaluationSites}{a 4-dimensional array containing the community evaluation metrics for the evaluation sites of each run (\code{NA} means that the site was used for calibaration)}
  \item{PA.allSites}{a 4-dimensional array of the binary prediction for all sites and runs under the different thresholding appraoches.}
}
\references{
Scherrer, D., D'Amen, M., Mateo, M.R.G., Fernandes, R.F. & Guisan , A. (2018) How to best threshold and validate stacked species assemblages? Community optimisation might hold the answer. Methods in Ecology and Evolution, in review
}
\author{
Daniel Scherrer <daniel.j.a.scherrer@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ecospat.CCV.createDataSplitTable}}; \code{\link{ecospat.CCV.communityEvaluation.prob}}
}
\examples{
\donttest{
#Loading species occurence data and remove empty communities
data(ecospat.testData)
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
                                     
#Thresholding all the predictions and calculating the community evaluation metrics
myCCV.communityEvaluation.bin <- ecospat.CCV.communityEvaluation.bin(
      ccv.modeling.data = myCCV.Models, 
      thresholds = c('MAX.KAPPA', 'MAX.ROC','PS_SDM'),
      community.metrics= c('SR.deviation','Sorensen'),
      parallel = FALSE,
      cpus = 4)

}
}
