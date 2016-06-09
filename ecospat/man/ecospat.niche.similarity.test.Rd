\name{ecospat.niche.similarity.test}

\alias{ecospat.niche.similarity.test}

\title{Niche Similarity Test}

\description{Run a niche similarity test (see Warren et al 2008) based on two species occurrence density grids.}

\usage{ecospat.niche.similarity.test (z1, z2, rep, one.sided=T, rand.type)}

\arguments{
  \item{z1}{Species 1 occurrence density grid created by \code{ecospat.grid.clim}.}
  \item{z2}{Species 2 occurrence density grid created by \code{ecospat.grid.clim}.}
  \item{rep}{The number of replications to perform.}
  \item{one.sided}{If true, perform a one-sided test, otherwise a two-sided test.}
    \item{rand.type}{Type of randomization on the density grids (1 or 2).}
}

\details{Compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches (z2.sim) as available in the range of z2 (z2$Z). z2.sim has the same pattern as z2 but the center is randomly translatated in the availabe z2$Z space and weighted by z2$Z densities.
If rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.

If you encounter a problem during your analyses, please first read the FAQ section of "Niche overlap" in http://www.unil.ch/ecospat/home/menuguid/ecospat-resources/tools.html

}

\value{a list with $obs = observed overlaps, $sim = simulated overlaps, $p.D = p-value of the test on D, $p.I = p-value of the test on I.}

\references{Broennimann, O., M.C. Fitzpatrick, P.B. Pearman, B. Petitpierre, L. Pellissier, N.G. Yoccoz, W. Thuiller, M.J. Fortin, C. Randin, N.E. Zimmermann, C.H. Graham and A. Guisan. 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. \emph{Global Ecology and Biogeography}, \bold{21}, 481-497.

Warren, D.L., R.E. Glor and M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}, 2868-2883.
}

\author{Olivier Broennimann \email{olivier.broennimann@unil.ch}}

\seealso{\code{\link{ecospat.grid.clim.dyn}}, \code{\link{ecospat.niche.equivalency.test}}}
