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

\value{A list with two elements:

\itemize{
  \item \code{dyn} A matrix, where each cell represents a pixel in the environment space of the study area (i.e., the ordination). Pixel values represent the classification of that cell, 0: unoccupied in z1 and z2; -1: expansion (occupied only in z2); 1: unfilling (occupied only in z1); and 2: stability (or overlap, occupied in both z1 and z2).
  \item \code{dynamic.index.w}, a vector containing the index values for expansion, stability, and unfilling.
}

Expansion is calculated as the number of cells occupied only in z2 divided by the total number of cells occupied in z2 (i.e., including cells occupied in both z1 and z2). This number is weighted by the occupancy density estimate for cells in z2.

Stability is calculated as the number of cells occupied in z2 and z1, divided by the total number of cells occupied in z2. This is also weighted by the occupancy density estimate for cells in z2.

Unfilling is calculated as the number of cells occupied only in z1, divided by the total number of cells occupied in z1 (including those shared with z2). This value is weighted by the occupancy density for cells in z1.

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

\seealso{\code{\link{ecospat.grid.clim.dyn}}}
