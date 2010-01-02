\name{is.outlier}
\alias{is.outlier}      

\title{Efficiency outlier summary for objects of class 'modlist' or 'replist'}

\description{
After having applied function \code{\link{KOD}} on objects of class 'modlist' or 'replist', \code{is.outlier}
 returns a vector of logicals for each run if they are efficiency outliers or not.  
}

\usage{
is.outlier(object)
}

\arguments{
  \item{object}{a model of class 'modlist' or 'replist'.}
}

\value{
A vector of logicals with the run names.
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
ml <- modlist(reps)
res <- KOD(ml)
## which runs are outliers?
outl <- is.outlier(res)
outl
which(outl)    
}

\keyword{models}
\keyword{nonlinear}
