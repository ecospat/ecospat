\name{ecospat-package}
\alias{ecospat-package}
\alias{ecospat}

\docType{package}

\title{Spatial Ecology Miscellaneous Methods}
\description{
Collection of methods, utilities and data sets for the support of spatial ecology analyses with a focus on pre-, core and post- modelling analyses of species distribution, niche quantification and community assembly. Specifically, 

\bold{-Pre-modelling:}

Spatial autocorrelation --> \code{ecospat.mantel.correlogram};

Variable selection --> \code{ecospat.npred};

Climate Analalogy --> \code{ecospat.climan, ecospat.mess} and \code{ecospat.plot.mess};

Phylogenetic diversity measures --> \code{ecospat.calculate.pd};

Biotic Interactions --> \code{ecospat.co-occurrences} and \code{ecospat.Cscore};

Minimum Dispersal routes --> \code{ecospat.mdr};

Niche Quantification --> \code{ecospat.grid.clim.dyn, ecospat.niche.equivalency.test,} 

\code{ecospat.niche.similarity.test, ecospat.plot.niche, ecospat.plot.niche.dyn,} 

\code{ecospat.plot.contrib, ecospat.niche.overlap, ecospat.plot.overlap.test,} 

\code{ecospat.niche.dyn.index} and \code{ecospat.shift.centroids};

Data Preparation --> \code{ecospat.caleval, ecospat.cor.plot, ecospat.makeDataFrame,} 

\code{ecospat.occ.desaggregation, ecospat.rand.pseudoabsences, ecospat.rcls.grd,} 

\code{ecospat.recstrat_prop, ecospat.recstrat_regl} and \code{ecospat.sample.envar};

\bold{-Core Niche Modelling:}

Model evaluation --> \code{ecospat.cv.glm, ecospat.permut.glm, ecospat.cv.gbm,} 

\code{ecospat.cv.me, ecospat.cv.rf, ecospat.boyce, ecospat.CommunityEval,} 

\code{ecospat.cohen.kappa, ecospat.max.kappa, ecospat.max.tss, ecospat.meva.table,} 

\code{ecospat.plot.kappa, ecospat.plot.tss} and \code{ ecospat.adj.D2.glm};

Spatial predictions and projections --> \code{ecospat.ESM.Modeling,} 

\code{ecospat.ESM.EnsembleModeling, ecospat.ESM.Projection, ecospat.ESM.EnsembleProjection,} 

\code{ecospat.SESAM.prr, ecospat.migclim, ecospat.binary.model, ecospat.Epred}

 and \code{ecospat.mpa};

Variable Importance --> \code{ecospat.maxentvarimport};

\bold{-Post Modelling:}

Variance Partition --> \code{ecospat.varpart};

Spatial predictions of species assemblages --> \code{ecospat.cons_Cscore};

Range size quantification --> \code{ecospat.rangesize} and 

\code{ecospat.occupied.patch};

The \code{ecospat} package was written by current and former members and collaborators of the ecospat group of Antoine Guisan, Department of Ecology and Evolution (DEE) & Institute of Earth Surface Dynamics (IDYST), University of Lausanne, Switzerland. 

}
\details{
\tabular{ll}{
Package: \tab ecospat\cr
Type: \tab Package\cr
Version: \tab 2.2.0\cr
Date: \tab 2017-11-22\cr
License: \tab GPL \cr
}

}
\author{Olivier Broennimann [aut],
             Valeria Di Cola [cre, aut],
             Blaise Petitpierre [ctb],
             Frank Breiner [ctb],
             Manuela D`Amen [ctb],
             Christophe Randin [ctb],
             Robin Engler [ctb],
             Wim Hordijk [ctb],
             Julien Pottier [ctb],
             Mirko Di Febbraro [ctb],
             Loic Pellissier [ctb],
             Dorothea Pio [ctb],
             Ruben Garcia Mateo [ctb],
             Anne Dubuis [ctb],
             Daniel Scherrer [ctb],
             Luigi Maiorano [ctb],
             Achilleas Psomas [ctb],
             Charlotte Ndiribe [ctb]
             Nicolas Salamin [ctb],
             Niklaus Zimmermann [ctb],
             Antoine Guisan [aut]
             }

\keyword{ package }
