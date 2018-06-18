\name{ecospat.meva.table}
\alias{ecospat.meva.table}

\title{
Model Evaluation For A Given Threshold Value
}
\description{
Calculates values of a series of different evaluations metrics for a model and for a given threshold value

}
\usage{
    ecospat.meva.table(Pred, Sp.occ, th)
}
\arguments{
    \item{Pred}{
        A vector of predicted probabilities
    }
    \item{Sp.occ}{
        A vector of binary observations of the species occurrence
    }
    \item{th}{
        Threshold used to cut the probability to binary values
    }
}

\value{
    A contingency table of observations and predicted probabilities of presence values, and a list of evaluation metrics for the selected threshold.
}
\author{
Antoine Guisan \email{antoine.guisan@unil.ch} with contributions of Luigi Maiorano \email{luigi.maiorano@gmail.com}
}
\references{
Pearce, J. and S. Ferrier. 2000. Evaluating the predictive performance of habitat models developed using logistic regression. \emph{Ecol. Model.}, \bold{133}, 225-245.
}

\seealso{
   \code{\link[ecospat]{ecospat.max.kappa}}, \code{\link[ecospat]{ecospat.max.tss}}, \code{\link[ecospat]{ecospat.plot.tss}}, \code{\link[ecospat]{ecospat.cohen.kappa}}, \code{\link[ecospat]{ecospat.plot.kappa}}

}
\examples{

Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris

meva <- ecospat.meva.table (Pred, Sp.occ, 0.39)
}
\keyword{file}
