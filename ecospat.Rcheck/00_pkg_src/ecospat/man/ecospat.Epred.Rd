\name{ecospat.Epred}

\alias{ecospat.Epred}

\title{Prediction Mean}

\description{Calculate the mean (or weighted mean) of several predictions.}

\usage{ecospat.Epred (x, w=rep(1,ncol(x)), th=0)}

\arguments{
  \item{x}{A dataframe object with SDM predictions.}
  \item{w}{Weight of the model, e.g. AUC. The default is 1.}
  \item{th}{Threshold used to binarize.}
}

\details{The Weighted Average consensus method utilizes pre-evaluation of the predictive performance of the single-models. In this approach, half (i.e. four) of the eight single-models with highest accuracy are selected first, and then a WA is calculated based on the pre-evaluated AUC of the single-models}

\value{A weighted mean binary transformation of the models.}

\author{Blaise Petitpierre \email{bpetitpierre@gmail.com}}

\references{
        
Boyce, M.S.,  P.R. Vernier, S.E. Nielsen and  F.K.A. Schmiegelow. 2002. Evaluating resource selection functions. \emph{Ecol. Model.}, \bold{157}, 281-300.

Marmion, M., M. Parviainen, M. Luoto, R.K. Heikkinen andW.  Thuiller. 2009. Evaluation of consensus methods in predictive species distribution modelling.  \emph{Diversity and Distributions}, \bold{15}, 59-69.

}

\examples{
x <- ecospat.testData[c(92,96)]
mean <- ecospat.Epred (x, w=rep(1,ncol(x)), th=0.5)
}