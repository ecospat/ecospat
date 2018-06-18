\name{ecospat.shift.centroids}

\alias{ecospat.shift.centroids}

\title{Draw Centroid Arrows}

\description{Draw arrows linking the centroid of the native and exotic (non-native) distribution (continuous line) and between native and invaded extent (dashed line).
}

\usage{ecospat.shift.centroids(sp1, sp2, clim1, clim2,col)}

\arguments{
  \item{sp1}{The scores of the species native distribution along the the two first axes of the PCA.}
  \item{sp2}{The scores of the species invasive distribution along the the two first axes of the PCA.}
  \item{clim1}{The scores of the entire native extent along the the two first axes of the PCA.}
  \item{clim2}{The scores of the entire invaded extent along the the two first axes of the PCA.}
  \item{col}{Colour of the arrow.}
}

\details{Allows to visualize the shift of the niche centroids of the species and the centroids of the climatic conditions in the study area. To compare invasive species niche, the arrow links the centroid of the native and inasive distribution (continuous line) and between native and invaded extent (dashed line). }

\value{Arrow on the overlap test plot}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

