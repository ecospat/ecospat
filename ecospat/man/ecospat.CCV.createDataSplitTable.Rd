\name{ecospat.CCV.createDataSplitTable}
\alias{ecospat.CCV.createDataSplitTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a DataSplitTable for usage in ecospat.ccv.modeling.
%%  ~~function to do ... ~~
}
\description{Creates a DataSplitTable with calibration and evaluation data either for cross-validation or repeated split sampling at the community level (i.e., across all species).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ecospat.CCV.createDataSplitTable(NbRunEval, 
                                 DataSplit,
                                 validation.method,
                                 NbSites,
                                 sp.data=NULL,
                                 minNbPresences=NULL,
                                 minNbAbsences=NULL,
                                 maxNbTry=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NbRunEval}{number of cross-validation or split sample runs}
  \item{DataSplit}{proportion (\%) of sites used for model calibration}
  \item{validation.method}{the type of \code{DataSplitTable} that should be created. Must be either \code{'cross-validation'} or \code{'split-sample'}}
    \item{NbSites}{number of total sites available. Is ignored if sp.data is provided.}
    \item{sp.data}{a data.frame where the rows are sites and the columns are species (values 1,0)}
    \item{minNbPresences}{the desired minimum number of Presences required in each run}
    \item{minNbAbsences}{the desired minimum number of Absences required in each run}
    \item{maxNbTry}{number of random tries allowed to create a fitting DataSplitTable}
}
\details{
If a \code{sp.data} data.frame with species presences and absences is provided the function tries to create a \code{DataSplitTable} which ensures that the maximum possible number of species can be modelled (according to the specified minimum presences and absences.)
}
\value{
  \item{DataSplitTable}{a matrix with \code{TRUE/FALSE} for each model run (\code{TRUE}=Calibration point, \code{FALSE}=Evaluation point)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Daniel Scherrer <daniel.j.a.scherrer@gmail.com>
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ecospat.CCV.modeling}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
#Creating a DataSplitTable for 200 sites, 25 runs with an 
#80/20 calibration/evaluation cross-validation

DataSplitTable <- ecospat.CCV.createDataSplitTable(NbSites = 200, 
                                                   NbRunEval=25, 
                                                   DataSplit=80, 
                                                   validation.method='cross-validation')
                                                   
#Loading species occurence data and remove empty communities
testData <- ecospat.testData[,c(24,34,43,45,48,53,55:58,60:63,65:66,68:71)]
sp.data <- testData[which(rowSums(testData)>0), sort(colnames(testData))]

#Creating a DataSplitTable based on species data directly
DataSplitTable <- ecospat.CCV.createDataSplitTable(NbRunEval = 20,
                                                   DataSplit = 70,
                                                   validation.method = "cross-validation",
                                                   NbSites = NULL,
                                                   sp.data = sp.data, 
                                                   minNbPresence = 15, 
                                                   minNbAbsences = 15, 
                                                   maxNbTry = 250)
}

}
