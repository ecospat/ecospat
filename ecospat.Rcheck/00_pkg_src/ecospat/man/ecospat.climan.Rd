\name{ecospat.climan}

\alias{ecospat.climan}

\title{A climate analogy setection tool for the modeling of species distributions}

\description{Assess climate analogy between a projection extent (p) and a reference extent (ref, used in general as the background to calibrate SDMs)}

\usage{ecospat.climan (ref, p)}

\arguments{
  \item{ref}{A dataframe with the value of the variables (i.e columns) for each point of the reference exent.}
  \item{p}{A dataframe with the value of the variables (i.e columns) for each point of the projection exent.}

}

\value{
Returns a vector. Values below 0 are novel conditions at the univariate level (similar to the MESS), values between 0 and 1 are analog and values above 1 are novel covariate condtions. 
For more information see Mesgeran et al. (2014)
}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

\references{
Mesgaran, M.B., R.D. Cousens and B.L. Webber. 2014. Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. \emph{Diversity & Distributions}, \bold{20}, 1147-1159.
}

\examples{
x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] #A reference dataset
ecospat.climan(ref,p)

}
