\name{ecospat.ESM.VarContrib}
\alias{ecospat.ESM.VarContrib}

\title{
Variable contribution in ESM
}
\description{
calculates the variable contribution of each variable and method in an ESM model
}
\usage{
    ecospat.ESM.VarContrib(ESM.modeling.output,
                           ESM_EF.output)
}

\arguments{
    \item{ESM.modeling.output}{\code{BIOMOD.formated.data} object returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
    \item{ESM_EF.output}{\code{BIOMOD.formated.data} object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
}
\details{
Calculates the ratio between sum of weights of bivariate models where a focal variable was used and sum of weights of bivariate models where the focal variable was not used. The ratio is corrected for the number of models with or without the focal variable. This ratio gives an indication on the proportional contribution of the variable in the final ensemble model. A value of higher than 1 indicate that the focal variable has a higher contribution than average.
In the case of multiple methods (e.g., GLM, GAM...), the contributions are counted per method. For ensemble model, the contributions are then weighted means (based on the weighting score as chosen in ecospat.ESM.EnsembleModeling()) of single methods
    
}

\value{
Returns a dataframe with contribution values (i.e. proportional contribution) by variable and model

}
\author{
Olivier Broennimann <Olivier.Broennimann@unil.ch> with contributions of Heidi Mod \email{Heidi.mod1@gmail.com} and Daniel Scherrer \email{daniel.j.a.scherrer@gmail.com}
}

\seealso{
\code{\link[ecospat]{ecospat.ESM.Modeling}}
}