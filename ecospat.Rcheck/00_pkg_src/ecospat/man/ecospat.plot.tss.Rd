\name{ecospat.plot.tss}
\alias{ecospat.plot.tss}

\title{
    Plot True skill statistic (TSS)
}
\description{
    Plots the values for True skill statistic (TSS) along different thresholds.
}
\usage{
    ecospat.plot.tss(Pred, Sp.occ)
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
 A plot of the TSS values along different thresholds.
}
\author{
Luigi Maiorano \email{luigi.maiorano@gmail.com}
}
\references{

Liu, C., P.M. Berry, T.P. Dawson, and R.G. Pearson. 2005. Selecting thresholds of occurrence in the prediction of species distributions. \emph{Ecography}, \bold{28}, 385-393.

Liu, C., M. White and G. Newell. 2013. Selecting thresholds for the prediction of species occurrence with presence-only data. \emph{Journal of Biogeography}, \emph{40}, 778-789.

}

\seealso{
   \code{\link[ecospat]{ecospat.meva.table}}, \code{\link[ecospat]{ecospat.max.tss}}, \code{\link[ecospat]{ecospat.plot.kappa}}, \code{\link[ecospat]{ecospat.cohen.kappa}}, \code{\link[ecospat]{ecospat.max.kappa}}
}
\examples{
Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
ecospat.plot.tss(Pred, Sp.occ)
}

\keyword{file}
