\name{repform}
\alias{repform}

\title{Formatting of PCR data with or without replicates for analysis}

\description{
  This function formats the PCR data prior to the model analysis. The columns of the PCR data, 
  that need to be in a dataframe or matrix, are defined by a factor vector in respect to 
  cycle data and fluorescence data of (possible) replicates.
}

\usage{
  repform(x, group)
}

\arguments{
  \item{x}{a dataframe or matrix containing the PCR data.}
  \item{group}{a vector defining the data structure in respect to PCR cycles and replicated data}
 }

\details{
  The function calls the \code{stack} function. The column containing the cycle number must be termed '0'.
}

\value{
  A dataframe that is formatted for use in the \code{multdrc} analysis.
}

\author{
Andrej-Nikolai Spiess & Christian Ritz
}

\examples{

data(reps)
### take the first five columns of the dataframe, 
### first column contains the cycle numbers,
###  the next four columns are replicates
repData <- repform(reps[,1:5], c(0,1,1,1,1))
m2 <- multdrc(values ~ Cycles, curve = Curve, data = repData, fct = l5())
pcrplot(m2, type = "errbar")

### take the first five columns of the dataframe, 
### first column contains the cycle numbers,
###  the next four columns are two pairs of replicates
repData <- repform(reps[,1:5], c(0,1,1,2,2))
m2 <- multdrc(values ~ Cycles, curve = Curve, data = repData, fct = l5())
pcrplot(m2, type = "all", col = 1:2)

}

\keyword{models}
\keyword{nonlinear}
