\name{ecospat.mantel.correlogram}

\alias{ecospat.mantel.correlogram}

\title{Mantel Correlogram}

\description{Investigate spatial autocorrelation of environmental covariables within a set of occurrences as a function of distance.}

\usage{ecospat.mantel.correlogram (dfvar, colxy, n, colvar, max, nclass, nperm)}

\arguments{
  \item{dfvar}{A dataframe object with the environmental variables.}
  \item{colxy}{The range of columns for x and y in df.}
  \item{n}{The number of random occurrences used for the test.}
  \item{colvar}{The range of columns for variables in df.}
  \item{max}{The maximum distance to be computed in the correlogram.}
  \item{nclass}{The number of classes of distances to be computed in the correlogram.}
  \item{nperm}{The number of permutations in the randomization process.}
}

\details{Requires ecodist library. Note that computation time increase tremendously when using more than 500 occurrences (n>500)}

\value{
Draws a plot with distance vs. the mantel r  value. Black circles indicate that the values are significative different from zero. White circles indicate non significant autocorrelation. The selected distance is at the first white circle where values are non significative different from cero. 
}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}


\references{
Legendre, P. and M.J. Fortin. 1989. Spatial pattern and ecological analysis. \emph{Vegetatio}, \bold{80}, 107-138.}

\seealso{\code{\link{mgram}}}

\examples{
ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, colvar=3:7, 
max=1000, nclass=10, nperm=100)
}