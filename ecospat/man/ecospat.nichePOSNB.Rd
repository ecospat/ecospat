\name{ecospat.nichePOSNB}
\alias{ecospat.nichePOSNB}

\title{
Niche position and niche breadth)
}
\description{
    ecospat.nichePOSNB calculates the niche breadth and niche position of taxa along environmental gradients from abundance data. 
}
\usage{
ecospat.nichePOSNB (df,colvar,colfreq)
}

\arguments{
    \item{df}{dataframe with (relative) abundance data. Each row correspond to an abundance.}
    \item{colvar}{The column(s) in df corresponding to environmental axe(s).}
    \item{colfreq}{The columns in df corresponding to taxa frequencies.}
 }
 
\details{
The function calculates niche position and niche breadth of taxa along one or multiple environmental axes. Niche position is calculated as the mean of the variable, weighted by the relative abundance of the species. Niche breadth is calculated as the standard deviation of each variable, weighted by the relative abundance of the species at each sampling site. If more than one variable is investigated, the values of niche position and niche breadth are averaged to obtain mean niche properties.
}

\value{
The function returns a matrix containing the average niche position and niche breadth of each taxa along each environmental axi. Additionally, if multiple environmental axes are considered, the mean niche and mean position of each taxa along these axes is returned by the function.
}

\author{Lucie Malard & Olivier Broennimann \email{lucie.malard@unil.ch, olivier.broennimann@unil.ch}}

\references{L.A. Malard, H.K. Mod, N. Guex, O. Broennimann, E. Yashiro, E. Lara, E.D.A. Mitchell, H. Niculita-Hirzel & A. Guisan. The ecological niche of soil bacterial, archaeal, fungal and protist communities along environmental gradients in the Alps. 2021. submitted.}

\examples{

data(ecospat.testPOSNB)
df<-ecospat.testPOSNB
ecospat.nichePOSNB(df,colvar=c(2),colfreq = 6:17) # 1 axes
ecospat.nichePOSNB(df,colvar=c(2:3),colfreq = 6:17) # 2 axes
ecospat.nichePOSNB(df,colvar=c(2:5),colfreq = 6:17) # 4 axes #
}
