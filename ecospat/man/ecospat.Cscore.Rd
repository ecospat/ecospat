\name{ecospat.Cscore}
\alias{ecospat.Cscore}

\title{
Pairwise co-occurrence Analysis with calculation of the C-score index.
}
\description{
The function tests for nonrandom patterns of species co-occurrence in a presence-absence matrix. It calculates the C-score index for the whole community and for each species pair. Null communities have column sum fixed.
}
\usage{
ecospat.Cscore (data, nperm, outpath, verbose = FALSE)
}
\arguments{
    \item{data}{A presence-absence dataframe for each species (columns) in each location or grid cell (rows). Column names (species names) and row names (sampling plots).}
    \item{nperm}{The number of permutation in the null model.}
    \item{outpath}{Path to specify where to save the results.}
    \item{verbose}{Boolean indicating whether to print progress output during calculation. Default is FALSE.}
}
\details{
This function allows to apply a pairwise null model analysis (Gotelli and Ulrich 2010) to a presence-absence community matrix to determine which species associations are significant across the study area. The strength of associations is quantified by the C-score index (Stone and Roberts 1990) and a 'fixed-equiprobable' null model algorithm is applied.
The format required for input databases: a plots (rows) x species (columns) matrix. Input matrices should have column names (species names) and row names (sampling plots). 
NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002). 

}
\value{
The function returns the C-score index for the observed community (ObsCscoreTot), p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the working directory where the same metrics are calculated for each species pair (only the table with species pairs with significant p-values is saved in this version)
}
\author{
Christophe Randin \email{	christophe.randin@wsl.ch} and Manuela D'Amen <manuela.damen@msn.com>
}
\references{
Gotelli, N.J. and D.J. McCabe. 2002. Species co-occurrence: a meta-analysis of JM Diamond's
assembly rules model. \emph{Ecology}, \bold{83}, 2091-2096.

Gotelli, N.J. and W. Ulrich. 2010. The empirical Bayes approach as a tool to identify non-random species associations. \emph{Oecologia}, \bold{162}, 463-477

Stone, L. and A. Roberts, A. 1990. The checkerboard score and species distributions. \emph{Oecologia}, \bold{85}, 74-79

}

\seealso{
   \code{\link{ecospat.co_occurrences}} and \code{\link{ecospat.cons_Cscore}}
}
\examples{
\dontrun{
data<- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
nperm <- 10000
outpath <- getwd()
Cscore<-ecospat.Cscore(data, nperm, outpath)

}}
