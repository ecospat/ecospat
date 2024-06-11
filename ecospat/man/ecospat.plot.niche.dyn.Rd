\name{ecospat.plot.niche.dyn}

\alias{ecospat.plot.niche.dyn}

\title{Plot of categories of Niche dynamics and Species Density}

\description{Plot niche categories of niche dynamics between two species densities created by \code{ecospat.grid.clim.dyn}.}

\usage{ecospat.plot.niche.dyn(z1, z2, intersection = 0, title = "", name.axis1 =
                   "Axis 1", name.axis2 = "Axis 2", interest = 1, col.abn
                   = "lightgreen", col.unf = "green", col.exp = "red",
                   col.stab = "blue", col.pio = "pink", col.NA = "grey",
                   colZ1 = "green3", colZ2 = "red3", transparency = 70,...)}

\arguments{
  \item{z1}{A gridclim object for the native distribution.}
  \item{z2}{A gridclim object for the invaded range.}
  \item{intersection}{The quantile of the environmental density used to delimit marginal climates. 
  If \code{intersection=0}, the analysis is performed at the intersection between native and invaded range.
  If \code{intersection=0.05}, the analysis is performed at the intersection of the 5th quantile of both native and invaded environmental densities}
  \item{title}{The title of the plot.}
  \item{name.axis1}{A label for the first axis.}
  \item{name.axis2}{A label for the second axis}
  \item{interest}{Choose which density to plot: if \code{interest=1}, plot native density, if \code{interest=2}, plot invasive density.}
  \item{col.abn}{The color used to depict the abandonment niche.}
  \item{col.unf}{The color used to depict the unfilling niche.}
  \item{col.exp}{The color used to depict the expansion niche.}
  \item{col.stab}{The color used to depict the stability niche.}
  \item{col.pio}{The color used to depict the pioneering niche.}
  \item{col.NA}{The color used to depict the environments outside of both niches.}
  \item{colZ1}{The color used to delimit the native extent.}
  \item{colZ2}{The color used to delimit the invaded extent.}
  \item{transparency}{A value between 0 and 100 to set the transparency level of the niche categories}
  \item{\dots}{Other graphical parameters.}
}

\details{
Using the default colors, the plot will show the niche stability in blue, niche expansion in red, and niche unfilling in green. 
The solid contour line indicates the extent of environmental conditions that exists in the native and invaded ranges; 
the dotted contour line indicates the quantile indicated by the \code{quant} argument.
The densities of occurrences are displayed using gray shading. This shading shows occurences in the native or
invaded range only, as determined by the value of the \code{interest} argument. 
}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}
