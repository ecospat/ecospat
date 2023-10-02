\name{ecospat.binary.model}

\alias{ecospat.binary.model}

\title{Generate Binary Models}

\description{Generate a binary map from a continuous model prediction.}

\usage{ecospat.binary.model (Pred, Threshold)}

\arguments{
  \item{Pred}{a \code{SpatRaster} object predicted suitabilities from a SDM prediction.}
  \item{Threshold}{A threshold to convert continous maps into binary maps (e.g. the output of the function \code{ecospat.mpa}() or use the \code{optimal.thresholds} from PresenceAbsence R package.}
}

\details{This function generates a binary model prediction (presence/absence) from an original model applying a threshold. The threshold could be arbitrary, or be based on the maximum acceptable error of false negatives (i.e. percentage of the presence predicted as absences, omission error).}

\value{The binary model prediction (presence/absence).}

\references{Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation}, \bold{24}: 38-49.

Engler, R., A Guisan and L. Rechsteiner. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. \emph{Journal of Applied Ecology}, \bold{41}, 263-274.

Liu, C., Berry, P. M., Dawson, T. P. and R. G. Pearson. 2005. Selecting thresholds of occurrence in the prediction of species distributions. \emph{Ecography}, \bold{28}, 385-393.

Jimenez-Valverde, A. and J.M.Lobo. 2007. Threshold criteria for conversion of probability of species presence to either-or presence-absence. \emph{Acta oecologica}, \bold{31}, 361-369.

Liu, C., White, M. and G. Newell. 2013. Selecting thresholds for the prediction of species occurrence with presence-only data. \emph{J. Biogeogr.}, \bold{40}, 778-789.

Freeman, E.A. and G.G. Moisen. 2008. A comparison of the performance of threshold criteria for binary classification in terms of predicted prevalence and kappa. \emph{Ecological Modelling}, \bold{217}, 48-58.
}

\seealso{
   \code{\link{ecospat.mpa}}, \code{\link[PresenceAbsence]{optimal.thresholds}}
}

\author{Ruben G. Mateo \email{rubeng.mateo@gmail.com} with contributions of Frank Breiner \email{frank.breiner@wsl.ch} and Flavien Collart}

\examples{

\donttest{
library(terra)
# coordinates of the plots
xy <- ecospat.testData[,2:3]

# environmental data
predictors <- terra::rast(system.file("extdata","ecospat.testEnv.tif",package="ecospat"))

env <- terra::extract(predictors,xy,ID=FALSE)  
spData <- cbind.data.frame(occ=ecospat.testData$Veronica_alpina,env)
mod <- glm(occ~ddeg0+I(ddeg0^2)+srad68+I(srad68^2),data=spData,family = binomial())

# predict to entire dataset
pred <- terra::predict(predictors,mod,type="response") 
  
### make binary maps
# use MPA to convert suitability to binary map

mpa.cutoff <- ecospat.mpa(pred,xy[spData$occ==1,],perc = 0.9) # 90% presences encompassed
pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)

plot(pred.bin.mpa)
points(xy[spData$occ==1,])
}
}
