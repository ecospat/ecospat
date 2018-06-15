\name{ecospat.rand.pseudoabsences}

\alias{ecospat.rand.pseudoabsences}

\title{Sample Pseudo-Absences}

\description{Randomly sample pseudoabsences from an environmental dataframe covering the study area.}

\usage{ecospat.rand.pseudoabsences (nbabsences, glob, colxyglob, colvar="all", 
presence, colxypresence, mindist)}

\arguments{
  \item{nbabsences}{The number of pseudoabsences desired.}
  \item{glob}{A two-column dataframe (or a vector) of the environmental values (in column) for background pixels of the whole study area (in row).}
  \item{colxyglob}{The range of columns for x and y in glob.}
  \item{colvar}{The range of columns for the environmental variables in glob. colvar="all" keeps all the variables in glob in the final dataframe. colvar=NULL keeps only x and y.}
  \item{presence}{A presence-absence dataframe for each species (columns) in each location or grid cell (rows).}
  \item{colxypresence}{The range of columns for x and y in presence.}
  \item{mindist}{The minimum distance from presences within wich pseudoabsences should not be drawn (buffer distance around presences).}
}

\value{A dataframe of random absences.}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\examples{
glob <- ecospat.testData[2:8]
presence <- ecospat.testData[c(2:3,9)]
presence <- presence[presence[,3]==1,1:2]
ecospat.rand.pseudoabsences (nbabsences=10, glob=glob, colxyglob=1:2, colvar = "all", 
presence= presence, colxypresence=1:2, mindist=20)
}