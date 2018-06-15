\name{ecospat.cons_Cscore}
\alias{ecospat.cons_Cscore}
\title{Constrained Co-Occurrence Analysis.}
\description{
Co-occurrence Analysis & Environmentally Constrained Null Models. The function tests for non-random patterns of species co-occurrence in a presence-absence matrix. It calculates the C-score index for the whole community and for each species pair. An environmental constraint is applied during the generation of the null communities.
}
\usage{ecospat.cons_Cscore(presence,pred,nperm,outpath)}
\arguments{
    \item{presence}{A presence-absence dataframe for each species (columns) in each location or grid cell (rows) Column names (species names) and row names (sampling plots).}
    \item{pred}{A dataframe object with SDM predictions. Column names (species names SDM) and row names (sampling plots).}
    \item{nperm}{The number of permutation in the null model.}
    \item{outpath}{Path to specify where to save the results.}
        }
\details{An environmentally constrained approach to null models will provide a more robust evaluation of species associations by facilitating the distinction between mutually exclusive processes that may shape species distributions and community assembly.
The format required for input databases: a plots (rows) x species (columns) matrix. Input matrices should have column names (species names) and row names (sampling plots).
NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002 - Ecology)}

\value{Returns the C-score index for the observed community (ObsCscoreTot), the mean of C-score for the simulated communities (SimCscoreTot), p.value (PValTot) and standardized effect size (SES.Tot). It also saves a table in the specified path where the same 
metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version). }

\author{Anne Dubuis \email{anne.dubuis@gmail.com} and Manuela D`Amen \email{manuela.damen@unil.ch}}

\references{Gotelli, N.J. and D.J. McCabe. 2002. Species co-occurrence: a meta-analysis of JM Diamond`s assembly rules model. \emph{Ecology}, \bold{83}, 2091-2096.

        Peres-Neto, P.R., J.D. Olden and D.A. Jackson. 2001. Environmentally constrained null models: site suitability as occupancy criterion. \emph{Oikos}, \bold{93}, 110-120.
        }

\examples{
\dontrun{
presence <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred <- ecospat.testData[c(73:92)]
nperm <- 10000
outpath <- getwd()
ecospat.cons_Cscore(presence, pred, nperm, outpath)
}}
