\name{KOD}
\alias{KOD}

\title{(K)inetic (O)utlier (D)etection according to Bar et al. (2003) or using Partitioning Around Medoids}

\description{
Identifies and/or removes qPCR runs whose efficiency differs significantly from all other runs of the same group.
Uses either outlier identification based on a normal distribution by iterating all replicates using the remaining samples as a training set
 (Leave-One-Out analysis) or PAM (Partitioning Around Medoids) as a different approach based on clustering. See 'Details'.
}

\usage{
KOD(object, method = c("bar", "pam"), 
    efftype = c("sliwin", "sigfit", "expfit"), train = TRUE, 
    remove = FALSE, alpha = 0.05, verbose = TRUE, ...)
}

\arguments{
  \item{object}{an object of class 'modlist' or 'replist'.}
  \item{method}{which method to use for outlier identification. Method from Bar et al. (2003) is default.}
  \item{efftype}{which method to use for calculation of qPCR efficiency. Either sliding window method, sigmoidal fitting or exponential fitting.}
  \item{train}{logical. If \code{TRUE}, sample is left out when calculating the p-values, if \code{FALSE}, sample is included. See Details.}
  \item{remove}{logical. If \code{TRUE}, outlier runs are removed and the object is updated. If \code{FALSE}, the individual qPCR runs are tagged as 'outliers' or not. See 'Details'.}
  \item{alpha}{The p-value threshold for identifying outliers by \code{method = "bar"}.}
  \item{verbose}{logical. If \code{TRUE}, all calculation steps and results are displayed on the console.} 
  \item{...}{any other parameters to be passed to \code{\link{pam}}, \code{\link{sliwin}}, \code{\link{efficiency}} or \code{\link{expfit}}.}
}

\details{
When using the KOD method according to Bar et al. (2003), outliers are defined by removing the sample from the replicate group and 
 testing it against the remaining training set by using a normal distribution:
\deqn{P^* = 2 * \lbrack1 - \Phi(\frac{e_i - \mu_{train}}{\sigma_{train}})\rbrack < 0.05}
In case of \code{method = "pam"}, outliers are identified by representing singletons when setting cluster size \code{k = 2} 
 in function \code{\link{pam}}. This is not a statistical approach but one based on the distance to the medoids, in which the outliers 
 exhibit a 'null' distance to the second medoid. Using this method might give completely other outliers as the default.
}

\value{
An object of the same class as in \code{object} that is 'tagged' in its name (**name**) if it is an outlier and also with an item \code{$outlier} with outlier information (see \code{\link{is.outlier}}). If \code{remove = TRUE}, the 
 outlier runs are removed (and the fitting updated in case of a 'replist').  
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
Function \code{\link{is.outlier}} to get an outlier summary.
}

\references{
Kinetic Outlier Detection (KOD) in real-time PCR.\cr
Bar T, Stahlberg A, Muszta A, Kubista M.\cr
\emph{Nucl Acid Res} (2003), \bold{31}, e105. 
}

\examples{
## on a 'modlist'
## F7.3 detected as outlier
ml <- modlist(reps, model = l5)
res <- KOD(ml)

## 'pam' detects no outliers
res2 <- KOD(ml, method = "pam") 

## on a 'replist',
## several outliers identified
rl <- replist(ml, group = gl(7, 4))
res2 <- KOD(rl)

## remove outliers and use
## plot matrix
res3 <- KOD(ml, remove = TRUE)
plot(res3, which = "single") 
}

\keyword{models}
\keyword{nonlinear}
