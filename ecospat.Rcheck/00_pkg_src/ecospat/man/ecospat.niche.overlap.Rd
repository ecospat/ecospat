\name{ecospat.niche.overlap}

\alias{ecospat.niche.overlap}

\title{Calculate Niche Overlap}

\description{Calculate the overlap metrics D and I based on two species occurrence density grids z1 and z2 created by \code{ecospat.grid.clim}.}

\usage{ecospat.niche.overlap (z1, z2, cor)}

\arguments{
  \item{z1}{Species 1 occurrence density grid created by \code{ecospat.grid.clim}.}
  \item{z2}{Species 2 occurrence density grid created by \code{ecospat.grid.clim}.}
  \item{cor}{Correct the occurrence densities of each species by the prevalence of the environments in their range (TRUE = yes, FALSE = no).}
}

\details{if cor=FALSE, the z$uncor objects created by \code{ecospat.grid.clim} are used to calculate the overlap. if cor=TRUE, the z$cor objects are used.

If you encounter a problem during your analyses, please first read the FAQ section of "Niche overlap" in http://www.unil.ch/ecospat/home/menuguid/ecospat-resources/tools.html
}

\value{Overlap values D and I. D is Schoener's overlap metric (Schoener 1970). I is a modified Hellinger metric(Warren et al. 2008)}

\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}, 481-497.

Schoener, T.W. 1968. Anolis lizards of Bimini: resource partitioning in a complex fauna. \emph{Ecology}, \bold{49}, 704-726.

Warren, D.L., R.E. Glor and M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}, 2868-2883.
}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\seealso{\code{\link{ecospat.grid.clim.dyn}}}
