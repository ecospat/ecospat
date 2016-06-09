\name{ecospat.max.tss}
\alias{ecospat.max.tss}

\title{
    Maximum TSS
}
\description{
    Calculates values for True skill statistic (TSS) along different thresholds, considering this time 0.01 increments (i.e. 99 thresholds).
}
\usage{
    ecospat.max.tss(Pred, Sp.occ)
}
\arguments{
    \item{Pred}{
        A vector of predicted probabilities
    }
    \item{Sp.occ}{
        A vector of binary observations of the species occurrence
    }
}

\value{
    Return values for TSS for 99 thresholds at 0.01 increments. 
}
\author{
Luigi Maiorano \email{luigi.maiorano@gmail.com} with contributions of Antoine Guisan \email{antoine.guisan@unil.ch}
}
\references{
Liu, C., P.M. Berry, T.P. Dawson, and R.G. Pearson. 2005. Selecting thresholds of occurrence in the prediction of species distributions. \emph{Ecography}, \bold{28}, 385-393.}

\seealso{
   \code{\link[ecospat]{ecospat.meva.table}}, \code{\link[ecospat]{ecospat.max.kappa}}, \code{\link[ecospat]{ecospat.plot.tss}}, \code{\link[ecospat]{ecospat.cohen.kappa}}, \code{\link[ecospat]{ecospat.plot.kappa}}
}
\examples{

Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
TSS100 <- ecospat.max.tss(Pred, Sp.occ)
}

\keyword{file}
