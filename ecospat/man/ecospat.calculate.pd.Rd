\name{ecospat.calculate.pd}

\alias{ecospat.calculate.pd}

\title{Calculate Phylogenetic Diversity Measures}

\description{Calculate all phylogenetic diversity measures listed in Schweiger et al., 2008 (see full reference below).}

\usage{ecospat.calculate.pd (tree, data, method="spanning", type="clade", root=FALSE, 
average=FALSE, verbose=FALSE)}

\arguments{
  \item{tree}{The phylogenetic tree}
  \item{data}{A presence or absence (binary) matrix for each species (columns) in each location or grid cell (rows)}
  \item{method}{The method to use. Options are "pairwise", "topology", and "spanning". Default is "spanning".}
  \item{type}{Phylogenetic measure from those listed in Schweiger et al 2008. Options are "Q", "P", "W", "clade", "species", "J", "F", "AvTD","TTD", "Dd". Default is "clade".}
  \item{root}{Phylogenetic diversity can either be rooted or unrooted. Details in Schweiger et al 2008. Default is FALSE.}
  \item{average}{Phylogenetic diversity can either be averaged or not averaged. Details in Schweiger et al 2008. Default is FALSE.}
  \item{verbose}{Boolean indicating whether to print progress output during calculation. Default is FALSE.}
}

\details{Given a phylogenetic tree and a presence/absence matrix this script calculates phylogenetic diversity of a group of species across a given set of grid cells or locations. The library "ape" is required to read the tree in R. Command is "read.tree" or "read.nexus".
Options of type:
"P" is a normalized mearure of "Q".
"clade" is "PDnode" when root= FALSE, and is "PDroot" ehn root =TRUE.
"species" is "AvPD".
}

\value{This function returns a list of phylogenetic diversity values for each of the grid cells in the presence/absence matrix}

\author{Nicolas Salamin \email{nicolas.salamin@unil.ch} and Dorothea Pio \email{Dorothea.Pio@fauna-flora.org}}

\references{
{Schweiger, O., S. Klotz, W. Durka and I. Kuhn. 2008. A comparative test of phylogenetic diversity indices. \emph{Oecologia}, \bold{157}, 485-495.}

{Pio, D.V., O. Broennimann, T.G. Barraclough, G. Reeves, A.G. Rebelo, W. Thuiller, A. Guisan and N. Salamin. 2011. Spatial predictions of phylogenetic diversity in conservation decision making. \emph{Conservation Biology}, \bold{25}, 1229-1239.}

{Pio, D.V., R. Engler, H.P. Linder, A. Monadjem, F.P.D. Cotterill, P.J. Taylor, M.C. Schoeman, B.W. Price, M.H. Villet, G. Eick, N. Salamin and A. Guisan. 2014. Climate change effects on animal and plant phylogenetic diversity in southern Africa. \emph{Global Change Biology}, \bold{20}, 1538-1549.}
}

\examples{

fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
library(ape)
tree <-ape::read.tree(fpath)
data <- ecospat.testData[9:52] 

pd <- ecospat.calculate.pd(tree, data, method = "spanning", type = "species", root = FALSE, 
average = FALSE, verbose = TRUE )

plot(pd)
}
