\name{ecospat.plot.niche.dyn}

\alias{ecospat.plot.niche.dyn}

\title{Niche Categories and Species Density}

\description{Plot niche categories and species density created by \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.plot.niche.dyn(z1, z2, quant, title = "",
               name.axis1 = "Axis 1", name.axis2 = "Axis 2",
               interest = 1, colz1 = "#00FF0050", colz2 = "#FF000050",
               colinter = "#0000FF50", colZ1 = "green3", colZ2 = "red3")}

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

\details{
Using the default colors, the plot will show the overlap area in blue,
niche expansion in red, and niche unfilling in green. A solid contour
line indicates the extent of habitat conditions that exist in the native
and invaded ranges; the dotted contour line indicates the quantile
indicated by the \code{quant} argument. i.e., if \code{quant} is set to
0.25, 75\% of available habitats are inside the dotted contour, and 25\%
of available habitats are between the dotted and solid contour lines.

The colored polygons show the distribution of occurrences as
presence-absence values. The densities of occurrences are displayed
using gray shading, with the darkest shading indicating the densest
grouping of occurrences. This shading shows occurences in the native or
invaded range only, as determined by the value of the \code{interest}
argument. 
}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}
