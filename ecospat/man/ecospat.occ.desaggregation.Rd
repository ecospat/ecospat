\name{ecospat.occ.desaggregation}

\alias{ecospat.occ.desaggregation}

\title{Species Occurrences Desaggregation}

\description{Remove species occurrences in a dataframe that are closer to each other than a specified distance threshold.}

\usage{ecospat.occ.desaggregation (dfvar, colxy, colvar=NULL, min.dist, plot=TRUE)}

\arguments{
  \item{xy}{A dataframe with xy-coordinates (x-column must be named 'x' and y-column 'y')}
  \item{min.dist}{The minimun distance between points in the sub-dataframe.}
  \item{by}{Grouping element in the dataframe (e.g. species, NULL)}
}

#\details{The number of occurrences kept and excluded is shown in a dynamic plot. At the end of the selection process, a plot showing the resulting distribution of occurrences will be plotted if argument plot=FALSE}

\value{A subset of df with the columns specified in colvar.}

\author {Frank Breiner} \email{frank.breiner@unil.ch} with contributions of {Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\examples{

\dontrun{
spp <- ecospat.testNiche
sp1 <- spp[1:32,1:2]

occ.sp1 <- ecospat.occ.desaggregation(xy=sp1, min.dist=500, by=NULL)
}
}
