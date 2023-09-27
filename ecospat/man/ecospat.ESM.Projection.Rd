\name{ecospat.ESM.Projection}
\alias{ecospat.ESM.Projection}

\title{
Ensemble of Small Models: Projects Simple Bivariate Models Into New Space Or Time
}
\description{
This function projects simple bivariate models on new.env
}
\usage{
    ecospat.ESM.Projection(ESM.modeling.output, 
                           new.env,
                           name.env,
                           parallel,
                           cleanup)
}
\arguments{
    \item{ESM.modeling.output}{
    \code{list} object returned by \code{\link{ecospat.ESM.Modeling}}
    }
    \item{new.env}{
A set of explanatory variables onto which models will be projected. It could be a \code{data.frame}, a \code{matrix}, or a \code{SpatRaster} object. Make sure the column names (\code{data.frame} or \code{matrix}) or layer Names (\code{SpatRaster}) perfectly match with the names of variables used to build the models in the previous steps.
    }
    \item{name.env}{
A name for the new.env object. If not specified (default) the name of the new.env object will be used. It is necessary to specify a unique name when projecting various new.env objects in a loop. 
    }
    \item{parallel}{
Logical. If TRUE, the parallel computing is enabled}
    \item{cleanup}{
\bold{No more available}. Please use terra::TmpFiles instead 
}
\details{
The basic idea of ensemble of small models (ESMs) is to model a species distribution based on small, simple models, for example all possible bivariate models (i.e. models that contain only two predictors at a time out of a larger set of predictors), and then combine all possible bivariate models into an ensemble (Lomba et al. 2010; Breiner et al. 2015).

The ESM set of functions could be used to build ESMs using simple bivariate models which are averaged using weights based on model performances (e.g. AUC) accoring to Breiner et al (2015). They provide full functionality of the approach described in Breiner et al. (2015).

The name of \code{new.env} must be a regular expression (see ?regex)
}
\value{
Returns the projections for all selected models (same as in \code{biomod2})
See \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} for details.

}
\author{ Frank Breiner \email{frank.breiner@wsl.ch} 

with contributions of Olivier Broennimann and Flavien Collart \email{olivier.broennimann@unil.ch}
}
\references{

Lomba, A., L. Pellissier, C.F. Randin, J. Vicente, F. Moreira, J. Honrado and A. Guisan. 2010. Overcoming the rare species modelling paradox: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation}, \bold{143},2647-2657.

Breiner F.T., A. Guisan, A. Bergamini and M.P. Nobis. 2015. Overcoming limitations of modelling rare species by using ensembles of small models. \emph{Methods in Ecology and Evolution}, \bold{6},1210-1218.

Breiner F.T., Nobis M.P., Bergamini A., Guisan A. 2018. Optimizing ensembles of small models for predicting the distribution of species with few occurrences. \emph{Methods in Ecology and Evolution}. \doi{10.1111/2041-210X.12957}
}

\seealso{
 \code{\link[ecospat]{ecospat.ESM.Modeling}}
}
