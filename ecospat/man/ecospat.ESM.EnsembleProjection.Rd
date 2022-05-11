\name{ecospat.ESM.EnsembleProjection}
\alias{ecospat.ESM.EnsembleProjection}

\title{
Ensemble of Small Models: Projects Calibrated ESMs Into New Space Or Time.
}
\description{
This function projects calibrated ESMs into new space or time.

}
\usage{
    ecospat.ESM.EnsembleProjection( ESM.prediction.output, 
                                    ESM.EnsembleModeling.output,
                                    chosen.models = 'all')
}
\arguments{
    \item{ESM.prediction.output}{a list object returned by \code{\link[ecospat]{ecospat.ESM.Projection}}
    }
    \item{ESM.EnsembleModeling.output}{
            a list object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
     \item{chosen.models}{a character vector (either 'all' or a sub-selection of model names, e.g. c(GLM, GBM)) to remove models from the ensemble. Default is 'all'.}
}

\details{
The basic idea of ensemble of small models (ESMs) is to model a species distribution based on small, simple models, for example all possible bivariate models (i.e. models that contain only two predictors at a time out of a larger set of predictors), and then combine all possible bivariate models into an ensemble (Lomba et al. 2010; Breiner et al. 2015).

The ESM set of functions could be used to build ESMs using simple bivariate models which are averaged using weights based on model performances (e.g. AUC) according to Breiner et al. (2015). They provide full functionality of the approach described in Breiner et al. (2015).
For projections only the full models (100% of the data) are used but not the split models which are solely used for model evaluation and model weighting.
    
For further details please refer to \code{\link[biomod2]{BIOMOD_EnsembleForecasting}}.

}
\value{
Returns the projections of ESMs for the selected single models and their ensemble (data frame or raster stack). ESM.projections \sQuote{projection files} are saved on the hard drive projection folder. This files are either an \code{array} or a \code{RasterStack} depending the original projections data type.
Load these created files to plot and work with them.

}
\author{
Frank Breiner \email{frank.breiner@wsl.ch}

}
\references{

Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.

Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.

Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
}

\seealso{
\code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}

}