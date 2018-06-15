\name{ecospat.plot.overlap.test}

\alias{ecospat.plot.overlap.test}

\title{Plot Overlap Test}

\description{Plot a histogram of observed and randomly simulated overlaps, with p-values of equivalency or similarity tests.}

\usage{ecospat.plot.overlap.test (x, type, title)}

\arguments{
  \item{x}{Object created by 
  
  \code{ecospat.niche.similarity.test} or \code{ecospat.niche.equivalency.test}.}
  \item{type}{Must be either \dQuote{D} or \dQuote{I}.}
  \item{title}{The title for the plot.}
}

\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}, 481-497.}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\seealso{\code{\link{ecospat.niche.similarity.test}}, \code{\link{ecospat.niche.equivalency.test}}}
