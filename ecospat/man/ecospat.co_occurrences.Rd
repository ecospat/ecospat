\name{ecospat.co_occurrences}

\alias{ecospat.co_occurrences}

\title{Species Co-Occurrences}

\description{Calculate an index of species co-occurrences.}

\usage{ecospat.co_occurrences (data)}

\arguments{
  \item{data}{A presence-absence matrix for each species (columns) in each location or grid cell (rows) or a matrix with predicted suitability values. }
}

\details{Computes an index of co-occurrences ranging from 0 (never co-occurring) to 1 (always co-occuring).}

\value{The species co-occurrence matrix and box-plot of the co-occurrence indices}


\author{Loic Pellissier \email{loic.pellissier@unifr.ch}}

\references{
Pellissier, L., K.A. Brathen, J. Pottier, C.F. Randin, P. Vittoz, A. Dubuis, N.G. Yoccoz, T. Alm, N.E. Zimmermann and A. Guisan. 2010. Species distribution models reveal apparent competitive and facilitative effects of a dominant species on the distribution of tundra plants. \emph{Ecography}, \bold{33}, 1004-1014.

Guisan, A. and N. Zimmermann. 2000. Predictive habitat distribution models in ecology. \emph{Ecological Modelling}, \bold{135}:147-186
}

\examples{
\dontrun{
matrix <- ecospat.testData[c(9:16,54:57)]
ecospat.co_occurrences (data=matrix)
}}