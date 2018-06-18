\name{ecospat.npred}

\alias{ecospat.npred}

\title{Number Of Predictors}

\description{Calculate the maximum number of predictors to include in the model with a desired correlation between predictors.}

\usage{ecospat.npred (x, th)}

\arguments{
  \item{x}{Correlation matrix of the predictors.}
  \item{th}{Desired threshold of correlation between predictors.}
}

\value{Returns the number of predictors to use.}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}


\examples{
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)
}