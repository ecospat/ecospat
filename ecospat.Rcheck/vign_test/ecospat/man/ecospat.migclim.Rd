\name{ecospat.migclim}

\alias{ecospat.migclim}

\title{Implementing Dispersal Into Species Distribution Models}

\description{Enables the implementation of species-specific dispersal constraints into projections of species distribution models under environmental change and/or landscape fragmentation scenarios.}

\usage{ecospat.migclim ()}

\details{The MigClim model is a cellular automaton originally designed to implement dispersal constraints into projections of species distributions under environmental change and landscape fragmentation scenarios.}

\references{Engler, R., W. Hordijk and A. Guisan. 2012. The MIGCLIM R package -- seamless integration of dispersal constraints into projections of species distribution models. \emph{Ecography}, \bold{35}, 872-878.

Engler, R. and A. Guisan. 2009. MIGCLIM: predicting plant distribution and dispersal in a changing climate. \emph{Diversity and Distributions}, \bold{15}, 590-601.

Engler, R., C.F. Randin, P. Vittoz, T. Czaka, M. Beniston, N.E. Zimmermann and A. Guisan. 2009. Predicting future distributions of mountain plants under climate change: does dispersal capacity matter? \emph{Ecography}, \bold{32}, 34-45.
}

\author{Robin Engler \email{robin.engler@gmail.com}, Wim Hordijk \email{wim@WorldWideWanderings.net} and Loic Pellissier \email{loic.pellissier@unifr.ch}}

\examples{
\dontrun{
ecospat.migclim()
### Some example data files can be downloaded from the following web page:
### http://www.unil.ch/ecospat/page89413.html
###
### Run the example as follows (set the current working directory to the
### folder where the example data files are located):
###
data(MigClim.testData)
### Run MigClim with a data frame type input.
n <- MigClim.migrate (iniDist=MigClim.testData[,1:3],
hsMap=MigClim.testData[,4:8], rcThreshold=500,
envChgSteps=5, dispSteps=5, dispKernel=c(1.0,0.4,0.16,0.06,0.03),
barrier=MigClim.testData[,9], barrierType="strong",
iniMatAge=1, propaguleProd=c(0.01,0.08,0.5,0.92),
lddFreq=0.1, lddMinDist=6, lddMaxDist=15,
simulName="MigClimTest", replicateNb=1, overWrite=TRUE,
testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)}

}