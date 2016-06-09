\name{ecospat.boyce}

\alias{ecospat.boyce}

\title{Calculate Boyce Index}

\description{Calculate the Boyce index as in Hirzel et al. (2006). The Boyce index is used to assess model performance.}

\usage{ecospat.boyce (fit, obs, nclass=0, window.w="default", res=100, PEplot=T)}

\arguments{
  \item{fit}{A vector containing the predicted suitability values}
  \item{obs}{A vector containing the predicted suitability values of the validation points (presence records)}
  \item{nclass}{The number of classes or vector with class thresholds. If \code{nclass=0}, the Boyce index is calculated with a moving window (see next parameters)}
  \item{window.w}{The width of the moving window (by default 1/10 of the suitability range)}
  \item{res}{The resolution of the moving window (by default 100 focals)}
  \item{PEplot}{If true, plot the predicted to expected ratio along the suitability class}
}

\details{The Boyce index only requires presences and measures how much model predictions differ from random distribution of the observed presences across the prediction gradients (Boyce et al. 2002). It is thus the most appropriate metric in the case of presence-only models. It is continuous and varies between -1 and +1. Positive values indicate a model which present predictions are consistent with the distribution of presences in the evaluation dataset, values close to zero mean that the model is not different from a random model, negative values indicate counter predictions, i.e., predicting poor quality areas where presences are more frequent (Hirzel et al. 2006).}

\value{Returns the predicted-to-expected ratio for each class-interval: F.ratio

Returns the Boyce index value: Spearman.cor

Creates a graphical plot of the predicted to expected ratio along the suitability class}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com} and Frank Breiner \email{frank.breiner@unil.ch}}

\references{
Boyce, M.S., P.R. Vernier, S.E. Nielsen and F.K.A. Schmiegelow. 2002. Evaluating resource selection functions. \emph{Ecol. Model.}, \bold{157}, 281-300.

Hirzel, A.H., G. Le Lay, V. Helfer, C. Randin and A. Guisan. 2006. Evaluating the ability of habitat suitability models to predict species presences. \emph{Ecol. Model.}, \bold{199}, 142-152.}

\examples{
obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.boyce (fit = ecospat.testData$glm_Saxifraga_oppositifolia , obs, nclass=0, 
window.w="default", res=100, PEplot=TRUE)
}