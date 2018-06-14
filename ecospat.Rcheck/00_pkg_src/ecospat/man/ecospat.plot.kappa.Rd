\name{ecospat.plot.kappa}
\alias{ecospat.plot.kappa}

\title{
Plot Kappa
}
\description{
Plots the values for Cohen's Kappa along different thresholds.
}
\usage{
    ecospat.plot.kappa(Pred, Sp.occ)
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
A plot of the Cohen's Kappa values along different thresholds.
}
\author{
Luigi Maiorano \email{luigi.maiorano@gmail.com} with contributions of Valeria Di Cola \email{valeria.dicola@unil.ch}.
}
\references{
Liu, C., P.M. Berry, T.P. Dawson, and R.G. Pearson. 2005. Selecting thresholds of occurrence in the prediction of species distributions. \emph{Ecography}, \bold{28}, 385-393.

Landis, J.R. and G.G. Koch. 1977. The measurement of observer agreement for categorical data. \emph{biometrics}, \bold{33},159-174.

}

\seealso{
   \code{\link[ecospat]{ecospat.meva.table}}, \code{\link[ecospat]{ecospat.max.tss}}, \code{\link[ecospat]{ecospat.plot.tss}}, \code{\link[ecospat]{ecospat.cohen.kappa}}, \code{\link[ecospat]{ecospat.max.kappa}}
}
\examples{


Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
ecospat.plot.kappa(Pred, Sp.occ)
}

\keyword{file}
