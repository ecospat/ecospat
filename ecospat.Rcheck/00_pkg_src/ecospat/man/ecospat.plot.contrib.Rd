\name{ecospat.plot.contrib}

\alias{ecospat.plot.contrib}

\title{Plot Variables Contribution}

\description{Plot the contribution of the initial variables to the analysis (i.e. correlation circle). Typically these are the eigen vectors and eigen values in ordinations.}

\usage{ecospat.plot.contrib (contrib, eigen)}

\arguments{
  \item{contrib}{A dataframe of the contribution of each original variable on each axis of the analysis, i.e. the eigen vectors. }
  \item{eigen}{A vector of the importance of the axes in the ordination, i.e. a vector of eigen values.}
}

\details{Requires ade4 library. If using \code{\link{princomp}} , use $loadings and $sdev of the princomp object. if using \code{\link{dudi.pca}}, use $li and $eig of the dudi.pca object.}

\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}, 481-497.}
 
\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\seealso{\code{\link{ecospat.plot.niche.dyn}},\code{\link{ecospat.plot.overlap.test}},
\code{\link{ecospat.niche.similarity.test}},\code{princomp}}
