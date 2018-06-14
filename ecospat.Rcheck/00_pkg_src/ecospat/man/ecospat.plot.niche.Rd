\name{ecospat.plot.niche}

\alias{ecospat.plot.niche}

\title{Plot Niche}

\description{Plot a niche z created by \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.plot.niche (z, title, name.axis1, name.axis2, cor=FALSE)}

\arguments{
  \item{z}{A gridclim object for the species distribution created by \code{ecospat.grid.clim.dyn}.}
  \item{title}{A title for the plot.}
  \item{name.axis1}{A label for the first axis.}
  \item{name.axis2}{A label for the second axis.}
  \item{cor}{Correct the occurrence densities of the species by the prevalence of the environments in its range (TRUE = yes, FALSE = no).}
}

\details{if z is bivariate, a bivariate plot of the niche of the species. if z is univariate, a histogram of the niche of the species}

\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}, 481-497.}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\seealso{\code{\link{ecospat.grid.clim.dyn}}}
