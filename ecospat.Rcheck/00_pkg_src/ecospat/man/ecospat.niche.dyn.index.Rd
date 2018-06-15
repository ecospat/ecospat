\name{ecospat.niche.dyn.index}

\alias{ecospat.niche.dyn.index}

\title{Niche Expansion, Stability, and Unfilling}

\description{Calculate niche expansion, stability and unfilling.}

\usage{ecospat.niche.dyn.index (z1, z2, intersection=NA)}

\arguments{
  \item{z1}{A gridclim object for the native distribution.}
  \item{z2}{A gridclim object for the invaded range.}
  \item{intersection}{The quantile of the environmental density used to remove marginal climates. If \code{intersection=NA}, the analysis is performed on the whole environmental extent (native and invaded). If \code{intersection=0}, the analysis is performed at the intersection between native and invaded range. If \code{intersection=0.05}, the analysis is performed at the intersection of the 5th quantile of both native and invaded environmental densities.}
}

\details{If you encounter a problem during your analyses, please first read the FAQ section of "Niche overlap" in http://www.unil.ch/ecospat/home/menuguid/ecospat-resources/tools.html}

\value{A list of dynamic indices: dynamic.index.w [expansion.index.w, stability.index.w, restriction.index.w]}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

\seealso{\code{\link{ecospat.grid.clim.dyn}}}
