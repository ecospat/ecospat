\name{ecospat.plot.niche.dyn}

\alias{ecospat.plot.niche.dyn}

\title{Niche Categories and Species Density}

\description{Plot niche categories and species density created by \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.plot.niche.dyn(z1, z2, quant=0.5, title = "",
               name.axis1 = "Axis 1", name.axis2 = "Axis 2",
               interest = 1, col.unf = "#00FF0050", col.exp = "#FF000050",
               col.stab = "#0000FF50", colZ1 = "green3", colZ2 = "red3")}

\arguments{
  \item{z1}{A gridclim object for the native distribution.}
  \item{z2}{A gridclim object for the invaded range.}
  \item{quant}{The quantile of the environmental density used to delimit marginal climates.}
  \item{title}{The title of the plot.}
  \item{name.axis1}{A label for the first axis.}
  \item{name.axis2}{A label for the second axis}
  \item{interest}{Choose which density to plot: if \code{interest=1}, plot native density, if \code{interest=2}, plot invasive density.}
  \item{col.unf}{The color used to depict unfilling area.}
  \item{col.exp}{The color used to depict expansion area.}
  \item{col.stab}{The color used to depict stability area.}
  \item{colZ1}{The color used to delimit the native extent.}
  \item{colZ2}{The color used to delimit the invaded extent.}

}

\details{
Using the default colors, the plot will show the niche stability in blue, niche expansion in red, and niche unfilling in green. 
The solid contour line indicates the extent of environmental conditions that exists in the native and invaded ranges; 
the dotted contour line indicates the quantile indicated by the \code{quant} argument.
The densities of occurrences are displayed using gray shading. This shading shows occurences in the native or
invaded range only, as determined by the value of the \code{interest} argument. 
}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}
