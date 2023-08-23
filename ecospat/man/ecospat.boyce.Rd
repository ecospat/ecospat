\name{ecospat.boyce}

\alias{ecospat.boyce}

\title{Calculate Boyce Index}

\description{Calculate the Boyce index as in Hirzel et al. (2006). The Boyce index is used to assess model performance.}

\usage{ecospat.boyce (fit, obs, nclass=0, window.w="default", res=100, PEplot = TRUE, 
rm.duplicate = TRUE, method = 'spearman')}

\arguments{
  \item{fit}{A vector, Raster-Layer or SpatRaster containing the predicted suitability values}
  \item{obs}{A vector containing the predicted suitability values or xy-coordinates (if "fit" is a Raster-Layer) of the validation points (presence records)}
  \item{nclass}{The number of classes or vector with class thresholds. If \code{nclass=0}, the Boyce index is calculated with a moving window (see next parameters)}
  \item{window.w}{The width of the moving window (by default 1/10 of the suitability range)}
  \item{res}{The resolution of the moving window (by default 100 focals)}
  \item{PEplot}{If \code{TRUE}, plot the predicted to expected ratio along the suitability class}
  \item{rm.duplicate}{If \code{TRUE}, only one of the successive duplicated P/E values are used to compute the correlation. See Details}
  \item{method}{Method used to compute the correlation. The default is \code{'spearman'} but \code{'pearson'} and \code{'kendall'} can also be used. See Details}


  rm.duplicate = TRUE, method = 'spearman'
}

\details{The Boyce index only requires presences and measures how much model predictions differ from random distribution of the observed presences across the prediction gradients (Boyce et al. 2002). It is thus an appropriate metric in the case of presence-only models. It is continuous and varies between -1 and +1. Positive values indicate a model which present predictions are consistent with the distribution of presences in the evaluation dataset, values close to zero mean that the model is not different from a random model, negative values indicate counter predictions, i.e., predicting poor quality areas where presences are more frequent (Hirzel et al. 2006).

Removing the successive duplicated P/E values (\code{rm.duplicate = TRUE}) focuses more on the discriminative aspect of the predictions, lowering the assessment of the evaluation of the model resolution (sensu Hirzel et al. 2006 p. 150). However, it seems that in the initial version, dupplicated values were not removed.

In the initial publication on the continuous Boyce index, the correlation was set to \code{method = 'spearman'}. However, using \code{method = 'kendall'} or \code{method = 'pearson'} might be more informative about the accuracy of the predictions.}

\value{The function returns a list that contains a vector F.ratio (the predicted-to-expected ratio for each class-interval) and a numeric Spearman.cor (the Boyce index value)}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com} and Frank Breiner \email{frank.breiner@unil.ch}}

\references{
Boyce, M.S., P.R. Vernier, S.E. Nielsen and F.K.A. Schmiegelow. 2002. Evaluating resource selection functions. \emph{Ecol. Model.}, \bold{157}, 281-300.

Hirzel, A.H., G. Le Lay, V. Helfer, C. Randin and A. Guisan. 2006. Evaluating the ability of habitat suitability models to predict species presences. \emph{Ecol. Model.}, \bold{199}, 142-152.}

\examples{
obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.boyce (fit = ecospat.testData$glm_Saxifraga_oppositifolia , obs, nclass=0, 
window.w="default", res=100, PEplot = TRUE)
}
