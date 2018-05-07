\name{ecospat.mess}

\alias{ecospat.mess}

\title{MESS}

\description{Calculate the MESS (i.e. extrapolation) as in Maxent.}

\usage{ecospat.mess (proj, cal, w="default")}

\arguments{
  \item{proj}{A dataframe object with x, y and environmental variables, used as projection dataset.}
  \item{cal}{A dataframe object with x, y and environmental variables, used as calibration dataset.}
  \item{w}{The weight for each predictor (e.g. variables importance in SDM).}
}

\details{Shows the variable that drives the multivariate environmental similarity surface (MESS) value in each grid cell.}

\value{
  \item{MESS}{The mess as calculated in Maxent, i.e. the minimal extrapolation values.}
  \item{MESSw}{The sum of negative MESS values corrected by the total number of predictors. If there are no negative values, MESSw is the mean MESS.}
  \item{MESSneg}{The number of predictors on which there is extrapolation.}
}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}. Modified by Daniel Scherrer \email{daniel.j.a.scherrer@gmail.com}}

\references{
Elith, J., M. Kearney and S. Phillips. 2010. The art of modelling range-shifting species. \emph{Methods in ecology and evolution}, \bold{1}, 330-342.
          }

\seealso{\code{\link{ecospat.plot.mess}}}


\examples{
x <- ecospat.testData[c(2,3,4:8)]
proj <- x[1:90,] #A projection dataset.
cal <- x[91:300,] #A calibration dataset

#Create a MESS object 
mess.object <- ecospat.mess (proj, cal, w="default")

#Plot MESS 
ecospat.plot.mess (mess.object, cex=1, pch=15)
}