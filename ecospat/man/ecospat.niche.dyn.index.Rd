\name{ecospat.niche.dyn.index}

\alias{ecospat.niche.dyn.index}

\title{Niche Pioneering, Expansion, Stability, Unfilling and Abandonment}

\description{Calculate niche Pioneering, Expansion, Stability, Unfilling and Abandonment.}

\usage{ecospat.niche.dyn.index (z1, z2, intersection=0)}

\arguments{
  \item{z1}{A gridclim object for the native distribution.}
  \item{z2}{A gridclim object for the invaded range.}
  \item{intersection}{
  The quantile of the environmental density used to remove marginal climates. 
  If \code{intersection=0}, the analysis is performed at the intersection between native and invaded range. 
  If \code{intersection=0.05}, the analysis is performed at the intersection of the 5th quantile of both native and invaded environmental densities.}
}

\details{
The function allows to perform niche dynamics analyses between two species. Originally based on Expansion, Stability and Unfilling indices as proposed in Guisan et al. 2014, it now includes the quantification of categories of niche dynamics proposed by Atwater et al. 2018, where the Pioneering niche represents the part of the niche Expansion in the invaded range that corresponds to environments not present in the native range and the Abandonment niche the part of the Unfilling niche in the native range that corresponds to environments not present in the invaded range. If you encounter a problem during your analyses, please first read the FAQ section of "Niche overlap" in https://www.unil.ch/ecospat/home/menuguid/ecospat-resources/tools.html
}

\value{A list of three objects: a rast object which classifies each pixels of the environmental space into Pioneering, Expansion, Stability, Unfilling or Abandonment Niche. A vector with the Expansion, Stability and Unfilling index values as proposed in Guisan et al. 2014. A vector with the count of pixels for each of the five niche dynamics categories.}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}, Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\seealso{\code{\link{ecospat.grid.clim.dyn}}}

\references{
Guisan, A., Petitpierre, B., Broennimann, O., Daehler, C., Kueffer, C. 2014. Unifying niche shift studies: insights from biological invasions. 
\emph{Trends in ecology and Ecolution}, \bold{29},260-269.

Atwater, D.Z., Ervine, C. and Barney J.N. 2018. Climatic niche shifts are common in introduced plants. \emph{Nature ecology and Evolution}, \bold{2}, 34-43.
}