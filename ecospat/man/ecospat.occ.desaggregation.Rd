\name{ecospat.occ.desaggregation}

\alias{ecospat.occ.desaggregation}

\title{Species Occurrences Desaggregation}

\description{Remove species occurrences in a dataframe that are closer to each other than a specified distance threshold.}

\usage{ecospat.occ.desaggregation (dfvar, colxy, colvar=NULL, min.dist, plot=TRUE)}

\arguments{
  \item{dfvar}{A dataframe with x, y, and variables.}
  \item{colxy}{The range of columns for x and y in df.}
  \item{colvar}{The range of columns for variables in df.}
  \item{min.dist}{The minimun distance threshold in the sub-dataframe.}
  \item{plot}{A boolean indicating whether to plot the resulting distribution of occurrences.}
}

\details{The number of occurrences kept and excluded is shown in a dynamic plot. At the end of the selection process, a plot showing the resulting distribution of occurrences will be plotted if argument plot=FALSE}

\value{A subset of df with the columns specified in colvar.}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\examples{

\dontrun{
spp <- ecospat.testNiche
sp1 <- spp[1:32,1:3]

occ.sp1 <- ecospat.occ.desaggregation(dfvar=sp1,colxy=2:3,colvar=NULL, min.dist=500,plot=TRUE)
}
}