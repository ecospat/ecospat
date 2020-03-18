\name{ecospat.mdr}
\alias{ecospat.mdr}

\title{
Minimum Dispersal Routes)
}
\description{
\code{ecospat.mdr} is a function that implement a minimum cost arborescence approach to analyse the invasion routes of a species from dates occurrence data.
}
\usage{
ecospat.mdr (data, xcol, ycol, datecol, mode, rep, mean.date.error, fixed.sources.rows)
}

\arguments{
    \item{data}{dataframe with occurence data. Each row correspond to an occurrence.}
    \item{xcol}{The column in data containing x coordinates.}
    \item{ycol}{The column in data containing y coordinates.}
    \item{datecol}{The column in data containing dates.}
    \item{mode}{"observed", "min" or "random". "observed" calculate routes using real dates. "min" reorder dates so the the total length of the routes are minimal. "random" reatribute dates randomly.}
    \item{rep}{number of iteration of the analyse. if > 1, boostrap support for each route is provided.}
    \item{mean.date.error}{mean number of years to substract to observed dates. It is the mean of the truncated negative exponential distribution from which the time to be substracted is randomly sampled.}
    \item{fixed.sources.rows}{the rows in data (as a vector) corresponding to source occurrence(s) that initiated the invasion(s). No incoming routes are not calculated for sources.}
}
\details{
The function draws an incoming route to each occurence from the closest occurrence already occupied (with a previous date) and allows to substract a random number of time (years) to the observed dates from a truncated negative exponential distribution. It is possible to run several iterations and to get boostrap support for each route.
\code{itexp} and \code{rtexp} functions are small internal functions to set a truncated negative exponential distribution.
}
\value{
A list is returned by the function with in positon [[1]], a datafame with each row corresponding to a route (with new/old coordinates, new/old dates, length of the route, timespan, dispersal rate), in position [[2]] the total route length, in position [[3]] the median dispersal rate and in position [[4]] the number of outgoing nodes (index of clustering of the network)
}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\references{
Hordijk, W. and O. Broennimann. 2012. Dispersal routes reconstruction and the minimum cost arborescence problem. \emph{Journal of theoretical biology}, \bold{308}, 115-122.

Broennimann, O., P. Mraz, B. Petitpierre, A. Guisan, and H. Muller-Scharer. 2014. Contrasting spatio-temporal climatic niche dynamics during the eastern and western invasions of spotted knapweed in North America.\emph{Journal of biogeography}, \bold{41}, 1126-1136.
}
\examples{
   \donttest{
library(maps)
data(ecospat.testMdr)
data<- ecospat.testMdr

fixed.sources.rows<-order(data$date)[1:2] #first introductions

#plot observed situation
plot(data[,2:1],pch=15,cex=0.5)
points(data[fixed.sources.rows,2:1],pch=19,col="red")
text(data[,2]+0.5,data[,1]+0.5,data[,3],cex=0.5)
map(add=T)

# mca 
obs<-ecospat.mdr(data=data,
xcol=2,
ycol=1,
datecol=3,
mode="min",
rep=100,
mean.date.error=1,
fixed.sources.rows)

#plot results
lwd<-(obs[[1]]$bootstrap.value)
x11();plot(obs[[1]][,3:4],type="n",xlab="longitude",ylab="latitude")
arrows(obs[[1]][,1],obs[[1]][,2],obs[[1]][,3],obs[[1]][,4],length = 0.05,lwd=lwd*2)
map(add=T)
points(data[fixed.sources.rows,2:1],pch=19,col="red")
text(data[fixed.sources.rows,2]+0.5,data[fixed.sources.rows,1]+0.5,data[fixed.sources.rows,3],
cex=1,col="red")
title(paste("total routes length : ",
round(obs[[2]],2)," Deg","\n","median dispersal rate : ",
round(obs[[3]],2)," Deg/year","\n","number of outcoming nodes : ",
obs[[4]]))
}
}
