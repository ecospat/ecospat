\name{ecospat.nicheNBmean}
\alias{ecospat.nicheNBmean}

\title{
Niche position and niche breadth)
}
\description{
    ecospat.nicheNBmean calculates the weighted mean niche breadth accross several axes from a nichePOSNB object. 
}
\usage{
ecospat.nichePOSNB (POSNB,w)
}

\arguments{
    \item{POSNB}{an object created by the function ecospat.nichePOSNB}
    \item{w}{a vector with a weight for each environmental axes.}
 }
 
\details{
The function calculates the weighted mean niche breadth. The weights w are rescaled so that their sum equals 1.
}

\value{
The function returns a 1 column dataframe with the weighted mean niche breadth for each taxa.
}

\author{Lucie Malard & Olivier Broennimann \email{lucie.malard@unil.ch, olivier.broennimann@unil.ch}}

\references{L.A. Malard, H.K. Mod, N. Guex, O. Broennimann, E. Yashiro, E. Lara, E.D.A. Mitchell, H. Niculita-Hirzel & A. Guisan. The ecological niche of soil bacterial, archaeal, fungal and protist communities along environmental gradients in the Alps. 2021. Accepted in Soil Biology and Biochemistry.}

\examples{

data(ecospat.testNichePOSNB)
df<-ecospat.testNichePOSNB

# 1 axes
POSNB<-ecospat.nichePOSNB(df,colvar=c(2),colfreq = 6:17) 

# 2 axes
POSNB<-ecospat.nichePOSNB(df,colvar=c(2:3),colfreq = 6:17) 
ecospat.nicheNBmean(POSNB,w=c(2,1))

# 4 axes
POSNB<-ecospat.nichePOSNB(df,colvar=c(2:5),colfreq = 6:17) 
ecospat.nicheNBmean(POSNB,w=c(1,0.8,0.2,0.1))

}
