\name{ecospat.plot.niche.dyn}

\alias{ecospat.plot.niche.dyn}

\title{Niche Categories and Species Density}

\description{Plot niche categories and species density created by \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.plot.niche.dyn (z1, z2, quant, title, 
name.axis1, name.axis2, interest, colz1, colz2,colinter, colZ1, colZ2)}

\arguments{
  \item{z1}{A gridclim object for the native distribution.}
  \item{z2}{A gridclim object for the invaded range.}
  \item{quant}{The quantile of the environmental density used to delimit marginal climates.}
  \item{title}{The title of the plot.}
  \item{name.axis1}{A label for the first axis.}
  \item{name.axis2}{A label for the second axis}
  \item{interest}{Choose which density to plot: if \code{interest=1}, plot native density, if \code{interest=2}, plot invasive density.}
  \item{colz1}{The color used to depict unfilling area.}
  \item{colz2}{The color used to depict expansion area.}
  \item{colinter}{The color used to depict overlap area.}
  \item{colZ1}{The color used to delimit the native extent.}
  \item{colZ2}{The color used to delimit the invaded extent.}

}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}
