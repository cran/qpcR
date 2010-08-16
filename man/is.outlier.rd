\name{is.outlier}
\alias{is.outlier}      

\title{Outlier summary for objects of class 'pcrfit, 'modlist' or 'replist'}

\description{
After using \code{\link{SOD}} or \code{\link{KOD}} on objects of class 'pcrfit', 'modlist' or 'replist', \code{is.outlier}
 returns a vector of logicals for each run if they are efficiency (\code{KOD}) or sigmoidal shape (\code{SOD}) outliers or not.  
}

\usage{
is.outlier(object)
}

\arguments{
  \item{object}{an object of class 'pcrfit', 'modlist' or 'replist' resulting from any of the outlier functions.}
}

\value{
A vector of logicals with run names.
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{KOD}} and \code{\link{SOD}}.
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

## test for sigmoidal outliers:
## create a non-amplification curve
## and throw into 'reps' dataset
dat <- reps
dat[, 10] <- rnorm(49, 1, 0.5)
ml <- modlist(dat)
res <- SOD(ml)
is.outlier(res)    
}

\keyword{models}
\keyword{nonlinear}
