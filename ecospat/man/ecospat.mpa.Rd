\name{ecospat.mpa}

\alias{ecospat.mpa}

\title{Minimal Predicted Area}

\description{Calculate the minimal predicted area.}

\usage{ecospat.mpa (Pred, Sp.occ.xy, perc)}

\arguments{
  \item{Pred}{Numeric or RasterLayer predicted suitabilities from a SDM prediction.}
  \item{Sp.occ.xy}{xy-coordinates of the species (if Pred is a RasterLayer).}
  \item{perc}{Percentage of Sp.occ.xy that should be encompassed by the binary map.}

}

\details{The minimal predicted area (MPA) is the minimal surface obtained by considering all pixels with predictions above a defined probability threshold (e.g. 0.7) that still encompasses 90 percent of the species` occurrences (Engler et al. 2004).}

\value{
  Returns the minimal predicted area.
}

\author{Frank Breiner \email{frank.breiner@wsl.ch}}

\references{
Engler, R., A. Guisan and L. Rechsteiner. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. \emph{Journal of Applied Ecology}, \bold{41}, 263-274.}


\examples{
data(ecospat.testData)
obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
[which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.mpa(obs)
ecospat.mpa(obs,perc=1) ## 100 percent of the presences encompassed
}