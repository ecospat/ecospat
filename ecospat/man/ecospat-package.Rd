\name{ecospat-package}
\alias{ecospat-package}
\alias{ecospat}

\docType{package}

\title{Spatial Ecology Miscellaneous Methods}

\description{
Collection of methods, utilities and data sets for the support of spatial ecology analyses with a focus on pre-, core and post- modelling analyses of species distribution, niche quantification and community assembly. The \code{ecospat} package was written by current and former members and collaborators of the ecospat group of Antoine Guisan, Department of Ecology and Evolution (DEE) & Institute of Earth Surface Dynamics (IDYST), University of Lausanne, Switzerland. 

\bold{Pre-modelling:}

\itemize{

\item Spatial autocorrelation:\code{\link{ecospat.mantel.correlogram}}

\item Variable selection: \code{\link{ecospat.npred}}

\item Climate Analalogy: \code{\link{ecospat.climan}, \link{ecospat.mess}, \link{ecospat.plot.mess}}

\item Phylogenetic diversity measures: \code{\link{ecospat.calculate.pd}}

\item Biotic Interactions: \code{\link{ecospat.cons_Cscore}, \link{ecospat.Cscore}, \link{ecospat.co_occurrences}} 

\item Minimum Dispersal routes: \code{\link{ecospat.mdr}}

\item Niche Quantification: \code{\link{ecospat.grid.clim.dyn}, \link{ecospat.niche.equivalency.test}, \link{ecospat.niche.similarity.test}, \link{ecospat.plot.niche}, \link{ecospat.plot.niche.dyn}, \link{ecospat.plot.contrib}, \link{ecospat.niche.overlap}, \link{ecospat.plot.overlap.test}, \link{ecospat.niche.dyn.index}, \link{ecospat.shift.centroids}, \link{ecospat.niche.dynIndexProjGeo}, \link{ecospat.niche.zProjGeo}}

\item Data Preparation: \code{\link{ecospat.caleval}, \link{ecospat.cor.plot}, \link{ecospat.makeDataFrame}, \link{ecospat.occ.desaggregation}, \link{ecospat.rand.pseudoabsences}, \link{ecospat.rcls.grd}, \link{ecospat.recstrat_prop}, \link{ecospat.recstrat_regl}, \link{ecospat.sample.envar}}
}

\bold{Core Niche Modelling:}

\itemize{

\item Model evaluation: \code{\link{ecospat.cv.glm}, \link{ecospat.permut.glm}, \link{ecospat.cv.gbm}, \link{ecospat.cv.me}, \link{ecospat.cv.rf}, \link{ecospat.boyce}, \link{ecospat.CommunityEval}, \link{ecospat.cohen.kappa}, \link{ecospat.max.kappa}, \link{ecospat.max.tss}, \link{ecospat.meva.table}, \link{ecospat.plot.kappa}, \link{ecospat.plot.tss}, \link{ecospat.adj.D2.glm}, \link{ecospat.CCV.createDataSplitTable}, \link{ecospat.CCV.modeling}, \link{ecospat.CCV.communityEvaluation.bin}, \link{ecospat.CCV.communityEvaluation.prob}}

\item Spatial predictions and projections: \code{\link{ecospat.ESM.Modeling}, \link{ecospat.ESM.EnsembleModeling}, \link{ecospat.ESM.Projection}, \link{ecospat.ESM.EnsembleProjection}, \link{ecospat.SESAM.prr}, \link{ecospat.binary.model}, \link{ecospat.Epred}, \link{ecospat.mpa}}

\item Variable Importance: \code{\link{ecospat.maxentvarimport}, \link{ecospat.ESM.VarContrib}}
}

\bold{Post Modelling:}

\itemize{
\item Variance Partition: \code{\link{ecospat.varpart}}

\item Spatial predictions of species assemblages: \code{\link{ecospat.cons_Cscore}}

\item Range size quantification: \code{\link{ecospat.rangesize}, \link{ecospat.occupied.patch}}
}
}

