\name{ecospat.cohen.kappa}
\alias{ecospat.cohen.kappa}

\title{Cohen's Kappa
}
\description{
Calculates Cohen's kappa and variance estimates, within a 95 percent confidence interval.
}
\usage{
    ecospat.cohen.kappa(xtab)
}
\arguments{
    \item{xtab}{
        A symmetric agreement table.
    }

}
\details{
The argument xtab is a contingency table.
xtab <- table(Pred >= th, Sp.occ)
}
\value{
A list with elements 'kap', 'vark', 'totn' and 'ci' is returned. 
'kap' is the cohen's kappa, 'vark' is the variance estimate within a 95 percent confidence interval, 'totn' is the number of plots and 'ci' is the confidence interval. 
}
\author{
Christophe Randin \email{christophe.randin@wsl.ch} with contributions of Niklaus. E. Zimmermann \email{niklaus.zimmermann@wsl.ch} and Valeria Di Cola \email{valeria.dicola@unil.ch}
}
\references{
Bishop, Y.M.M., S.E. Fienberg and P.W. Holland. 1975. Discrete multivariate analysis: Theory and Practice. Cambridge, MA: MIT Press. pp. 395-397.

Pearce, J. and S. Ferrier. 2000. Evaluating the predictive performance of habitat models developed using logistic regression. \emph{Ecol. Model.}, \bold{133}, 225-245.
}

\seealso{
   \code{\link[ecospat]{ecospat.meva.table}}, \code{\link[ecospat]{ecospat.max.tss}}, \code{\link[ecospat]{ecospat.plot.tss}}, \code{\link[ecospat]{ecospat.plot.kappa}}, \code{\link[ecospat]{ecospat.max.kappa}}
}
\examples{
Pred <- ecospat.testData$glm_Agrostis_capillaris
Sp.occ <- ecospat.testData$Agrostis_capillaris
th <- 0.39 # threshold
xtab <- table(Pred >= th, Sp.occ)

ecospat.cohen.kappa(xtab)
}
