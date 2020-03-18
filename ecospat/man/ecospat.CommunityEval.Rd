\name{ecospat.CommunityEval}

\alias{ecospat.CommunityEval}

\title{Community Evaluation}

\description{Calculate several indices of accuracy of community predictions.}

\usage{ecospat.CommunityEval (eval, pred, proba, ntir, verbose = FALSE)}

\arguments{
  \item{eval}{A matrix of observed presence-absence (ideally independent from the dataset used to fit species distribution models) of the species with n rows for the sites and s columns for the species.}
  \item{pred}{A matrix of predictions for the s species in the n sites. Should have the same dimension as eval.}
  \item{proba}{Logical variable indicating whether the prediction matrix contains presences-absences (FALSE) or probabilities (TRUE).}
  \item{ntir}{Number of trials of presence-absence predictions if pred is a probability matrix.}
  \item{verbose}{Boolean indicating whether to print progress output during calculation. Default is FALSE.}
}

\details{This function calculates several indices of accuracy of community predictions based on stacked predictions of species ditribution models. In case proba is set to FALSE the function returns one value per index and per site. In case proba is set to TRUE the function generates presences-absences based on the predicted probabilities and returns one value per index, per site and per trial.}

\value{A list of evaluation metrics calculated for each site (+ each trial if proba is set to TRUE):

  deviance.rich.pred: the deviation of the predicted species richness to the observed
  
  overprediction: the proportion of species predicted as present but not observed among the species predicted as present
    
  underprediction: the proportion of species predicted as absent but observed among the species observed as present
    
  prediction.success: the proportion of species correctly predicted as present or absent
  
  sensitivity: the proportion of species correctly predicted as present among the species observed as present
  
  specificity : the proportion of species correctly predicted as absent among the species observed as absent
  
  kappa: the proportion of specific agreement
  
  TSS: sensitivity+specificity-1
  
  similarity: the similarity of community composition between the observation and the prediction. The calculation is based on the Sorenses index.
  
  Jaccard: this index is a widely used metric of community similarity.}

\references{Pottier, J., A. Dubuis, L. Pellissier, L. Maiorano, L. Rossier, C.F. Randin, P. Vittoz and A. Guisan. 2013. The accuracy of plant assemblage prediction from species distribution models varies along environmental gradients. \emph{Global Ecology and Biogeography}, \bold{22}, 52-63.}

\author{Julien Pottier \email{julien.pottier@clermont.inra.fr} 

with contribution of Daniel Scherrer \email{daniel.scherrer@unil.ch}, Anne Dubuis \email{anne.dubuis@gmail.com} 
and Manuela D'Amen \email{manuela.damen@unil.ch}}

\examples{
\donttest{

data(ecospat.testData)
eval <- ecospat.testData[c(53,62,58,70,61,66,65,71,69,43,63,56,68,57,55,60,54,67,59,64)]
pred <- ecospat.testData[c(73:92)]

ecospat.CommunityEval (eval, pred, proba=TRUE, ntir=10)
}}