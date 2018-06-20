\name{ecospat.adj.D2.glm}
\alias{ecospat.adj.D2.glm}

\title{Calculate An Adjusted D2}

\description{This function is used for calculating an adjusted D2 from a calibrated GLM object}

\usage{
    ecospat.adj.D2.glm(glm.obj)
}
\arguments{
    \item{glm.obj}{Any calibrated GLM object with a binomial error distribution}}

\details{This function takes a calibrated GLM object with a binomial error distribution and returns an evaluation of the model fit.
The measure of the fit of the models is expressed as the percentage of explained deviance adjusted by the number of degrees
of freedom used (similar to the adjusted-R2 in the case of Least-Square regression; see Weisberg 1980) and is called the adjusted-D2
 (see guisan and Zimmermann 2000 for details on its calculation).
}

\value{Returns an adjusted D square value (proportion of deviance accounted for by the model).}

\author{
Christophe Randin \email{christophe.randin@unibas.ch} and Antoine Guisan \email{antoine.guisan@unil.ch}
}
\references{
Weisberg, S. 1980. Applied linear regression. Wiley.

Guisan, A., S.B. Weiss and A.D. Weiss. 1999. GLM versus CCA spatial modeling of plant species distribution. \emph{Plant Ecology}, \bold{143}, 107-122.

Guisan, A. and N.E. Zimmermann. 2000. Predictive habitat distribution models in ecology. \emph{Ecol. Model.}, \bold{135}, 147-186.
}

\examples{

data(ecospat.testData)
glm.obj<-glm(Achillea_millefolium~ddeg+mind+srad+slp+topo, 
family = binomial, data=ecospat.testData)

ecospat.adj.D2.glm(glm.obj)

}
