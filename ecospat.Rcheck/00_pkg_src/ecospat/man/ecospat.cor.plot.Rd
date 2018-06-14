\name{ecospat.cor.plot}
\alias{ecospat.cor.plot}

\title{Correlation Plot}

\description{A scatter plot of matrices, with bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal. Useful for descriptive statistics of small data sets (better with less than 10 variables).
}

\usage{ecospat.cor.plot(data)
}
\arguments{
    \item{data}{A dataframe object with environmental variables.}
}

\details{Adapted from the pairs help page. Uses panel.cor, and panel.hist, all taken from the help pages for pairs. It is a simplifies version of \code{ pairs.panels}() function of the package \code{psych}.
}
\value{A scatter plot matrix is drawn in the graphic window. The lower off diagonal draws scatter plots, the diagonal histograms, the upper off diagonal reports the Pearson correlation.}

\author{Adjusted by L. Mathys, 2006, modified by N.E. Zimmermann}

\examples{
data <- ecospat.testData[,4:8]
ecospat.cor.plot(data)
}
