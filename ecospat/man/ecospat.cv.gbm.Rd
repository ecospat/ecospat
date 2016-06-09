\name{ecospat.cv.gbm}

\alias{ecospat.cv.gbm}

\title{GBM Cross Validation}

\description{K-fold and leave-one-out cross validation for GBM.}

\usage{ecospat.cv.gbm (gbm.obj, data.cv, K=10, cv.lim=10, jack.knife=F)}

\arguments{
  \item{gbm.obj}{A calibrated GBM object with a binomial error distribution. Attention: users have to tune model input parameters according to their study!}
  \item{data.cv}{A dataframe object containing the calibration data set with the same names for response and predictor variables.}
  \item{K}{Number of folds. 10 is recommended; 5 for small data sets.}
  \item{cv.lim}{Minimum number of presences required to perform the K-fold cross-validation.}
  \item{jack.knife}{If TRUE, then the leave-one-out / jacknife cross-validation is performed instead of the 10-fold cross-validation.}
}

\details{This function takes a calibrated GBM object with a binomial error distribution and returns predictions from a stratified 10-fold cross-validation or a leave-one-out / jack-knived cross-validation. Stratified means that the original prevalence of the presences and absences in the full dataset is conserved in each fold.}

\value{Returns a dataframe with the observations (obs) and the corresponding predictions by cross-validation or jacknife.}

\references{Randin, C.F., T. Dirnbock, S. Dullinger, N.E. Zimmermann, M. Zappa and A. Guisan. 2006. Are niche-based species distribution models transferable in space? \emph{Journal of Biogeography}, \bold{33}, 1689-1703.
                     

Pearman, P.B., C.F. Randin, O. Broennimann, P. Vittoz, W.O. van der Knaap, R. Engler, G. Le Lay, N.E. Zimmermann and A. Guisan. 2008. Prediction of plant species distributions across six millennia. \emph{Ecology Letters}, \bold{11}, 357-369.}

\author{Christophe Randin \email{christophe.randin@unibas.ch} and Antoine Guisan \email{antoine.guisan@unil.ch}}


\examples{
\dontrun{
gbm <- ecospat.cv.gbm (gbm.obj= get ("gbm.Agrostis_capillaris", envir=ecospat.env), 
ecospat.testData, K=10, cv.lim=10, jack.knife=FALSE)
}
}