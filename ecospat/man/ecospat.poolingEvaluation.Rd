\name{ecospat.poolingEvaluation}
\encoding{UTF-8}
\alias{ecospat.poolingEvaluation}

\title{
  Evaluation of species distribution models using the pooling procedure
}
\description{
  This function evaluates species distribution models using 100\% of the dataset by pooling the different runs of the cross validation as in Collart et al. 2021
}
\usage{
  ecospat.poolingEvaluation(fit,
                            calib,
                            resp,
                            AlgoName = NULL,
                            metrics = c("SomersD","AUC","MaxTSS","MaxKappa","Boyce"),
                            ensembleEvaluation=FALSE,
                            w=NULL,
                            metricToEnsemble = "MaxTSS")

}
\arguments{
  \item{fit}{
    a list containing \emph{n} data.frame or matrix, where \emph{n} corresponds to the number of algorithm you want to evaluate. The data.frames (matrices) need to contain the model predictions (ranging between 0 and 1) resulting from the different runs of cross-validation. These data.frame need to have the same number of rows as in the full dataset (100\% of the occurrences and 100\% of the absences or background points) and a number of column equal to the number of cross-validation runs
}
  \item{calib}{
    a logical matrix with a number of rows equal to the full dataset and a number of column corresponding to the number of cross-validation runs. The value TRUE is to mention the elements that where used to calibrate the models whereas FALSE corresponds to the one that will be used for the evaluation (\emph{NB} the points used to calibrate the models during a cross-validation run should be the same across algoritms)
  }
  \item{resp}{
  a numeric vector where 1 corresponds to a species response and 0 to an absence (or background point) with a length corresponding to number of rows in calib
  }
  \item{AlgoName}{
  a character vector for giving a name to each elements of fit. If NULL, the position in the list will be used instead.
  }
  \item{metrics}{a vector of evaluation metrics chosen among "SomersD", "AUC", "MaxTSS", "MaxKappa", "Boyce"}
    \item{ensembleEvaluation}{
  logical. If TRUE, the ensemble model will be evaluated applying a weighted mean across algorithms.  
  }
  \item{w}{
  a numeric vector of the weights used to realize the ensemble model. The length should match the number of algorithms.
  }
    \item{metricToEnsemble}{
  character. Metric to use to ensemble the models with a weighted mean when w is not given. The metric should be one in metrics}
  }

\details{
Because a minimum sample size is needed to evaluate models (see Jiménez-Valverde, 2020), this function uses the approach from Collart et al.(2021), which consists to pool the suitability values of the hold-out data (evaluation dataset) across replicates. As the same data point (presence or absence or background point) is presumably sampled in several replicates, the suitability values for each data point is consequently averaged across replicates where they were sampled. This procedure generates a series of independent suitability values with a size approximately equal (as some data points may not have been sampled by chance in any of the \emph{n} replicates) to that of the number of data point.
}
\value{
  a list containing:
  \item{evaluations}{
  a matrix with the evaluation scores based on the different modelling algorithms and based on the consensus across the modelling algorithms (called here "ensemble")
  }
  \item{fit}{
  a matrix of predicted values resulting from the pooling procedure and used to compute the evaluation scores. The column \emph{resp} is where the species occurs or not
  }
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
  \code{\link[ecospat]{ecospat.ESM.EnsembleEvaluation}}
}
\examples{
  set.seed(42)
  resp <- c(rep(1,15),rep(0,85)) #15 presences and 85 absences
  #Generating a fake fit object whith two algorithms and 3 cross-vaidation
  fit <- list(matrix(0,nc=3,nr=100),
              matrix(0,nc=3,nr=100))
  fit[[1]][1:15,] = sample(seq(0,1, by=0.01),15*3,prob=c(rep(1,51),rep(10,50)),replace=TRUE)
  fit[[2]][1:15,] = sample(seq(0,1, by=0.01),15*3,prob=c(rep(1,51),rep(10,50)),replace=TRUE)
  fit[[1]][16:100,] = sample(seq(0,1, by=0.01),85*3,prob=c(rep(10,51),rep(1,50)),replace=TRUE)
  fit[[2]][16:100,] = sample(seq(0,1, by=0.01),85*3,prob=c(rep(10,51),rep(1,50)),replace=TRUE)
  
  # Generating a calib object where 80\% of the dataset is used to calibrate the model 
  # and 20\% to evaluate it
  calib <- matrix(TRUE,nc=3,nr=100)
  calib[c(sample(1:15,3),sample(16:100,17)),1]=FALSE 
  calib[c(sample(1:15,3),sample(16:100,17)),2]=FALSE 
  calib[c(sample(1:15,3),sample(16:100,17)),3]=FALSE 
  
  # Evaluation via the pooling procedure
  eval <- ecospat.poolingEvaluation(fit=fit,calib=calib,resp=resp,metrics=c("AUC","MaxTSS"))
  eval$evaluations
  
  # Evaluation including the ensemble model based on a weighted mean using MaxTSS
  evalEns <- ecospat.poolingEvaluation(fit=fit,calib=calib,resp=resp,ensembleEvaluation=TRUE,
                                       metrics=c("AUC","MaxTSS"))
  evalEns$evaluations
  
  # Evaluation including the ensemble model based on a mean by giving the same weight for 
  # each algorithm
  evalEns <- ecospat.poolingEvaluation(fit=fit,calib=calib,resp=resp,ensembleEvaluation=TRUE,
                                       metrics=c("AUC","MaxTSS"),w=c(1,1))
  evalEns$evaluations
}
