\name{ecospat.occ.desaggregation}

\alias{ecospat.occ.desaggregation}

\title{Species Occurrences Desaggregation}

\description{Remove species occurrences in a dataframe that are closer to each other than a specified distance threshold.}

\usage{ecospat.occ.desaggregation (xy, min.dist, by)}

\arguments{
  \item{xy}{A dataframe with xy-coordinates (x-column must be named 'x' and y-column 'y')}
  \item{min.dist}{The minimun distance between points in the sub-dataframe.}
  \item{by}{Grouping element in the dataframe (e.g. species, NULL)}
}

\details{This function will desaggregate the original number of occurrences, according to a specified distance.}

\value{A subset of the initial dataframe. Indicated the "initial", "kept" and "out number of points.}

\author{Frank Breiner \email{frank.breiner@unil.ch} 

with contributions of Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\examples{

\dontrun{
spp <- ecospat.testNiche
sp1 <- spp[1:32,2:3]

occ.sp1 <- ecospat.occ.desaggregation(xy=sp1, min.dist=500, by=NULL)
}
}