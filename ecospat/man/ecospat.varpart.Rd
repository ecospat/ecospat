\name{ecospat.varpart}

\alias{ecospat.varpart}

\title{Variation Partitioning For GLM Or GAM}

\description{Perform variance partitioning for binomial GLM or GAM based on the deviance of two groups or predicting variables.}

\usage{ecospat.varpart (model.1, model.2, model.12)}

\arguments{
  \item{model.1}{GLM / GAM calibrated on the first group of variables.}
  \item{model.2}{GLM / GAM calibrated on the second group of variables.}
  \item{model.12}{GLM / GAM calibrated on all variables from the two groups.}
}

\details{The deviance is calculated with the adjusted geometric mean squared improvement rescaled for a maximum of 1.}

\value{Return the four fractions of deviance as in Randin et al. 2009: partial deviance of model 1 and 2, joined deviance and unexplained deviance.}

\author{Christophe Randin \email{christophe.randin@unibas.ch}, Helene Jaccard and Nigel Gilles Yoccoz}

\references{Randin, C.F., H. Jaccard, P. Vittoz, N.G. Yoccoz and A. Guisan. 2009. Land use improves spatial predictions of mountain plant abundance but not presence-absence. \emph{Journal of Vegetation Science}, \bold{20}, 996-1008.}


\examples{
\donttest{
ecospat.cv.example()
ecospat.varpart (model.1= get ("glm.Achillea_atrata", envir=ecospat.env), 
model.2= get ("glm.Achillea_millefolium", envir=ecospat.env), 
model.12= get ("glm.Achillea_millefolium", envir=ecospat.env))
}
}