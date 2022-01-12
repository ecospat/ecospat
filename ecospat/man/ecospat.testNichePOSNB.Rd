\name{ecospat.testNichePOSNB}

\docType{data}

\alias{ecospat.testNichePOSNB}

\title{Test AVS Dataset For The Ecospat package}

\description{The dataset is contains frequencies of 15 bacterial consortium AVS for 16s site in the Western Swiss Alps along with 4 PCA scores representing environmental axes.}

\format{ecospat.testTree is a  16 rows (sites) x 18 collumns (14 AVS + 4 PCA axes) dataframe }

\references{L.A. Malard, H.K. Mod, N. Guex, O. Broennimann, E. Yashiro, E. Lara, E.D.A. Mitchell, H. Niculita-Hirzel & A. Guisan. The ecological niche of soil bacterial, archaeal, fungal and protist communities along environmental gradients in the Alps. 2021. Accepted in Soil Biology and Biochemistry.}

\author{Lucie Malard & Olivier Broennimann \email{lucie.malard@unil.ch, olivier.broennimann@unil.ch}}

\examples{
data(ecospat.testNichePOSNB)
df<-ecospat.testNichePOSNB
ecospat.nichePOSNB(df,colvar=c(2),colfreq = 6:17) # 1 axes
ecospat.nichePOSNB(df,colvar=c(2:3),colfreq = 6:17) # 2 axes
ecospat.nichePOSNB(df,colvar=c(2:5),colfreq = 6:17) # 4 axes #
}