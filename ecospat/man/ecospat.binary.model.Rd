\name{ecospat.binary.model}

\alias{ecospat.binary.model}

\title{Generate Binary Models}

\description{Generate a binary map from a continuous model prediction.}

\usage{ecospat.binary.model (Pred, Threshold)}

\arguments{
  \item{Pred}{RasterLayer predicted suitabilities from a SDM prediction.}
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

\author{Ruben G. Mateo \email{rubeng.mateo@gmail.com} with contributions of Frank Breiner \email{frank.breiner@wsl.ch}}

\examples{

\dontrun{
library(dismo)

### make a maxent model

# copy maxent.jar file in the right folder
path.from<-system.file("extdata", "maxent.jar", package="ecospat")
path.to <- paste0(system.file(package="dismo"), "/java/maxent.jar")
file.copy(path.from,path.to,overwrite = TRUE)

if (file.exists(path.to) & require(rJava)) {

  # get predictor variables
  fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
                     pattern='grd', full.names=TRUE )
  predictors <- stack(fnames)
  #plot(predictors)

  # file with presence points
  occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
  occ <- read.table(occurence, header=TRUE, sep=',')[,-1]
  colnames(occ) <- c("x","y")
  occ <- ecospat.occ.desaggregation(occ,min.dist=1)

  # fit a domain model, biome is a categorical variable
  me <- maxent(predictors, occ, factors='biome')

  # predict to entire dataset
  pred <- predict(me, predictors) 

  plot(pred)
  points(occ)

### convert suitability to binary map

mpa.cutoff <- ecospat.mpa(pred,occ,perc = 0.9) # 90 percent of the presences encompassed
pred.bin.mpa <- ecospat.binary.model(pred,mpa.cutoff)

plot(pred.bin.mpa)
points(occ)
}
}
}
