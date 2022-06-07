\name{ecospat.cv.glm}

\alias{ecospat.cv.glm}

\title{GLM Cross Validation}

\description{K-fold and leave-one-out cross validation for GLM.}

\usage{ecospat.cv.glm (glm.obj, K=10, cv.lim=10, jack.knife = FALSE, verbose = FALSE)}

\arguments{
  \item{glm.obj}{Any calibrated GLM object with a binomial error distribution.}
  \item{K}{Number of folds. 10 is recommended; 5 for small data sets.}
  \item{cv.lim}{Minimum number of presences required to perform the K-fold cross-validation.}
  \item{jack.knife}{If TRUE, then the leave-one-out / jacknife cross-validation is performed instead of the 10-fold cross-validation.}
  \item{verbose}{Boolean indicating whether to print progress output during calculation. Default is FALSE.}
}

\details{This function takes a calibrated GLM object with a binomial error distribution and returns predictions from a stratified 10-fold cross-validation or a leave-one-out / jack-knived cross-validation. Stratified means that the original prevalence of the presences and absences in the full dataset is conserved in each fold.}

\value{Returns a dataframe with the observations (obs) and the corresponding predictions by cross-validation or jacknife.}

\references{Randin, C.F., T. Dirnbock, S. Dullinger, N.E. Zimmermann, M. Zappa and A. Guisan. 2006. Are niche-based species distribution models transferable in space? \emph{Journal of Biogeography}, \bold{33}, 1689-1703.

Pearman, P.B., C.F. Randin, O. Broennimann, P. Vittoz, W.O. van der Knaap, R. Engler, G. Le Lay, N.E. Zimmermann and A. Guisan. 2008. Prediction of plant species distributions across six millennia. \emph{Ecology Letters}, \bold{11}, 357-369.}

\author{Christophe Randin \email{christophe.randin@unibas.ch} and Antoine Guisan \email{antoine.guisan@unil.ch}}


\examples{
\donttest{

if(require("rms",quietly=TRUE)){
  data('ecospat.testData')

  # data for Soldanella alpina
  data.Solalp<- ecospat.testData[c("Soldanella_alpina","ddeg","mind","srad","slp","topo")] 

  # gbm model for Soldanella alpina

  glm.Solalp <- glm(Soldanella_alpina ~ pol(ddeg,2) + pol(mind,2) + pol(srad,2) + pol(slp,2) 
      + pol(topo,2), data = data.Solalp, family = binomial)

  # cross-validated predictions
  glm.pred <- ecospat.cv.glm (glm.obj=glm.Solalp , K=10, cv.lim=10, jack.knife=FALSE)
}
}
}