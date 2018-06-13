\name{ecospat.max.kappa}
\alias{ecospat.max.kappa}

\title{
    Maximum Kappa
}
\description{
Calculates values for Cohen's Kappa along different thresholds, considering this time 0.01 increments (i.e. 99 thresholds).
}
\usage{
    ecospat.max.kappa(Pred, Sp.occ)
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
    Return values for Cohen's Kappa for 99 thresholds at 0.01 increments. 
}
\author{
Antoine Guisan \email{antoine.guisan@unil.ch} with contributions of Luigi Maiorano \email{luigi.maiorano@gmail.com} and Valeria Di Cola \email{valeria.dicola@unil.ch}.
}
\references{
Liu, C., P.M. Berry, T.P. Dawson, and R.G. Pearson. 2005. Selecting thresholds of occurrence in the prediction of species distributions. \emph{Ecography}, \bold{28}, 385-393.

}

\seealso{
   \code{\link[ecospat]{ecospat.meva.table}}, \code{\link[ecospat]{ecospat.max.tss}}, \code{\link[ecospat]{ecospat.plot.tss}}, \code{\link[ecospat]{ecospat.cohen.kappa}}, \code{\link[ecospat]{ecospat.plot.kappa}}
}
\examples{

   \dontrun{
Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
kappa100 <- ecospat.max.kappa(Pred, Sp.occ)
   }

}
