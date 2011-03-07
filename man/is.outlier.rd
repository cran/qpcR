\name{is.outlier}
\alias{is.outlier}      

\title{Outlier summary for objects of class 'modlist' or 'replist'}

\description{
For model lists of class 'modlist' or 'replist', \code{is.outlier} returns a vector of logicals for each run if they are outliers (i.e. sigmoidal or kinetic) or not.  
}

\usage{
is.outlier(object)
}

\arguments{
  \item{object}{an object of class 'modlist' or 'replist'.}
}

\value{
A vector of logicals with run names.
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{KOD}}.
}

\examples{
## analyze in respect to amplification
## efficiency outliers
ml1 <- modlist(reps, check = "uni2")
res1 <- KOD(ml1)

## which runs are outliers?
outl <- is.outlier(res1)
outl
which(outl)

## test for sigmoidal outliers
## with the 'testdat' dataset
ml2 <- modlist(testdat, check = "uni2")
is.outlier(ml2)    
}

\keyword{models}
\keyword{nonlinear}
