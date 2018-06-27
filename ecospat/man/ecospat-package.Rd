\name{ecospat-package}
\alias{ecospat-package}
\alias{ecospat}

\docType{package}

\title{Spatial Ecology Miscellaneous Methods}

\description{
Collection of methods, utilities and data sets for the support of spatial ecology analyses with a focus on pre-, core and post- modelling analyses of species distribution, niche quantification and community assembly. The \code{ecospat} package was written by current and former members and collaborators of the ecospat group of Antoine Guisan, Department of Ecology and Evolution (DEE) & Institute of Earth Surface Dynamics (IDYST), University of Lausanne, Switzerland. 

\bold{-Pre-modelling:}

Spatial autocorrelation:
  \code{\link{ecospat.mantel.correlogram}}

Variable selection --> \code{\link{ecospat.npred}}

Climate Analalogy --> \code{\link{ecospat.climan, ecospat.mess, ecospat.plot.mess}}

Phylogenetic diversity measures --> \code{\link{ecospat.calculate.pd}}

Biotic Interactions --> \code{\link{ecospat.co-occurrences,ecospat.Cscore}}

Minimum Dispersal routes --> \code{\link{ecospat.mdr}}

Niche Quantification --> \code{\link{ecospat.grid.clim.dyn, ecospat.niche.equivalency.test,ecospat.niche.similarity.test, ecospat.plot.niche, ecospat.plot.niche.dyn,ecospat.plot.contrib, ecospat.niche.overlap, ecospat.plot.overlap.test,ecospat.niche.dyn.index,ecospat.shift.centroids}}

Data Preparation --> \code{\link{ecospat.caleval, ecospat.cor.plot, ecospat.makeDataFrame,ecospat.occ.desaggregation, ecospat.rand.pseudoabsences, ecospat.rcls.grd,ecospat.recstrat_prop, ecospat.recstrat_regl,ecospat.sample.envar}}

\bold{-Core Niche Modelling:}}

Model evaluation --> \code{\link{ecospat.cv.glm, ecospat.permut.glm, ecospat.cv.gbm,ecospat.cv.me, ecospat.cv.rf, ecospat.boyce, ecospat.CommunityEval,ecospat.cohen.kappa, ecospat.max.kappa, ecospat.max.tss, ecospat.meva.table, ecospat.plot.kappa, ecospat.plot.tss, ecospat.adj.D2.glm}}

Spatial predictions and projections --> \code{\link{ecospat.ESM.Modeling,ecospat.ESM.EnsembleModeling, ecospat.ESM.Projection, ecospat.ESM.EnsembleProjection,ecospat.SESAM.prr, ecospat.migclim, ecospat.binary.model, ecospat.Epred, ecospat.mpa}}

Variable Importance --> \code{\link{ecospat.maxentvarimport}}

\bold{-Post Modelling:}}

Variance Partition --> \code{\link{ecospat.varpart}}

Spatial predictions of species assemblages --> \code{\link{ecospat.cons_Cscore}}

Range size quantification --> \code{\link{ecospat.rangesize, ecospat.occupied.patch}}

}

