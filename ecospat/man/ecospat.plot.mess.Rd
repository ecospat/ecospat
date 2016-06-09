\name{ecospat.plot.mess}

\alias{ecospat.plot.mess}

\title{Plot MESS}

\description{Plot the MESS extrapolation index onto the geographical space.}

\usage{ecospat.plot.mess (xy, mess.object, cex=1, pch=15)}

\arguments{
  \item{xy}{The x and y coordinates of the projection dataset.}
  \item{mess.object}{A dataframe as returned by the \code{ecospat.mess} function.}
  \item{cex}{Specify the size of the symbol.}
  \item{pch}{Specify the point symbols.}
}
\value{Returns a plot of the the MESS extrapolation index onto the geographical space.} 

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

\references{
Elith, J., M. Kearney and S. Phillips. 2010. The art of modelling range-shifting species. \emph{Methods in ecology and evolution}, \bold{1}, 330-342.}

\seealso{\code{\link{ecospat.mess}}}


\examples{
\dontrun{
x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset.
cal <- x[91:300,] #A calibration dataset

#Create a MESS object 
mess.object <- ecospat.mess (proj, cal, w="default")

#Plot MESS 
ecospat.plot.mess (xy=proj[c(1:2)], mess.object, cex=1, pch=15)
}
}