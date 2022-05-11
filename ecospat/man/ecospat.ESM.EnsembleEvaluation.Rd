\name{ecospat.ESM.EnsembleEvaluation}
\alias{ecospat.ESM.EnsembleEvaluation}

\title{
Ensemble of Small Models: Evaluation via the Pooling procedure
}
\description{
This function evaluates the Ensemble of Small Models by pooling the different runs of the cross validation as in Collart et al. 2021
}
\usage{
ecospat.ESM.EnsembleEvaluation(ESM.modeling.output,
                               ESM.EnsembleModeling.output,
                               metrics = c("SomersD","AUC","MaxTSS","MaxKappa","Boyce"), 
                               EachSmallModels = F)
}
\arguments{
  \item{ESM.modeling.output}{a \code{list} returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
  \item{ESM.EnsembleModeling.output}{a \code{list} returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
  \item{metrics}{a vector of evaluation metrics chosen among "SomersD", "AUC", "MaxTSS", "MaxKappa", "Boyce"}
  \item{EachSmallModels}{should the individual bivariate models be evaluated by the pooling procedure?}
}
\details{
Because a minimum sample size is needed to evaluate models (see Jiménez-Valverde, 2020), this function uses the approach from Collart et al.(2021), which consists to pool the suitability values of the hold-out data (evaluation dataset) across replicates. As the same data point (presence or absence or background point) is presumably sampled in several replicates, the suitability values for each data point is consequently averaged across replicates where they were sampled. This procedure generates a series of independent suitability values with a size approximately equal (as some data points may not have been sampled by chance in any of the \emph{n} replicates) to that of the number of data point.
}
\value{
a list containing:
  \item{ESM.evaluations}{a matrix with the evaluation scores for the ESMs based on the different modelling algorithms and based on the consensus across the modelling algorithms (called here "ensemble")}
  \item{ESM.fit}{a matrix of predicted values resulting from the pooling procedure and used to compute the evaluation scores. The column \emph{resp} is where the species occurs or not}
  \item{ESM.evaluations.bivariate.models}{a list containing a matrix of evaluation scores for each bivariate models (generated only if EachSmallModels = T)}
  \item{ESM.fit.bivariate.models}{a list containing a matrix of of predicted values resulting from the pooling procedure for each bivariate models (generated only if EachSmallModels = T)}
}
\references{
Collart, F., Hedenäs, L., Broennimann, O., Guisan, A. and Vanderpoorten, A. 2021. Intraspecific differentiation: Implications for niche and distribution modelling. \emph{Journal of Biogeography}. \bold{48}, 415-426. \doi{10.1111/jbi.14009}

Jiménez-Valverde, A. 2020. Sample size for the evaluation of presence-absence models. \emph{Ecological Indicators}. \bold{114}, 106289. \doi{10.1016/j.ecolind.2020.106289}
}
\author{
Flavien Collart \email{flavien.collart@unil.ch}

with contributions of Olivier Broennimann \email{olivier.broennimann@unil.ch}
}

\seealso{
\code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}
}