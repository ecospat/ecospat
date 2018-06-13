\name{ecospat.caleval}

\alias{ecospat.caleval}

\title{Calibration And Evaluation Dataset}

\description{Generate an evaluation and calibration dataset with a desired ratio of disaggregation.}

\usage{ecospat.caleval (data, xy, row.num=1:nrow(data), nrep=1, ratio=0.7, 
disaggregate=0, pseudoabs=0, npres=0, replace=FALSE)}

\arguments{
  \item{data}{A vector with presence-absence (0-1) data for one species.}
  \item{xy}{The x and y coordinates of the projection dataset.}
  \item{row.num}{Row original number}
  \item{nrep}{Number of repetitions}
  \item{ratio}{Ratio of disaggregation}
  \item{disaggregate}{Minimum distance of disaggregation (has to be in the same scale as xy)}
  \item{pseudoabs}{Number of pseudoabsences}
  \item{npres}{To select a smaller number of presences from the dataset to be subsetted. The maximum number is the total number of presences}
  \item{replace}{F to replace de pseudoabsences}
}

\details{This functions generates two list, one with the calibration or training dataset and other list with the evaluation or testing dataset disaggregated with a minimum distance.}

\value{list("eval"=eval,"cal"=cal))}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}


\examples{
data <- ecospat.testData
caleval <- ecospat.caleval (data = ecospat.testData[53], xy = data[2:3], row.num = 1:nrow(data), 
nrep = 2, ratio = 0.7, disaggregate = 0.2, pseudoabs = 100, npres = 10, replace = FALSE)
caleval
}