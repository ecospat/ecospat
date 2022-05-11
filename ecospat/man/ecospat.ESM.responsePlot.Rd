\name{ecospat.ESM.responsePlot}
\alias{ecospat.ESM.responsePlot}

\title{
Produce response plots for ESMs
}
\description{
This function creates response plots (evaluation strips) for Ensmebles of Small Models (ESMs).

}
\usage{
    ecospat.ESM.responsePlot( ESM.EnsembleModeling.output,
                              ESM.modeling.output,
                              fixed.var.metric = 'median')
}
\arguments{
    \item{ESM.modeling.output}{a list object returned by \code{\link[ecospat]{ecospat.ESM.Modeling}}}
    \item{ESM.EnsembleModeling.output}{a list object returned by \code{\link[ecospat]{ecospat.ESM.EnsembleModeling}}}
    \item{fixed.var.metric}{either 'median' (default), 'mean', 'min' or 'max' specifying the statistic used to fix as constant the remaining variables when the predicted response is estimated for one of the variables. (same as in \code{\link[biomod2]{response.plot2}}} 
}

\details{
This function plots the response curves of a model for each variable, while keeping the remianing variables constant. This is an adaptation of the Evaluation Strip method proposed by Elith et al.(2005)
}

\value{
A plot of the response curves is produced (red line Ensemble, grey lines single algorithms) and a list with the output is provided.
}

\author{
Frank Breiner \email{frank.breiner@unil.ch}

}
\references{
Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting predicted responses from species distribution models. Ecological Modelling 186, 280-289.
  }

\seealso{
\code{\link[ecospat]{ecospat.ESM.Modeling}}
}