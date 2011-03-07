\name{SOD}
\alias{SOD}

\title{(S)igmoidal (O)utlier (D)etection using first/second derivative maxima}

\description{
Identifies and/or removes qPCR runs that lack sigmoidal structure (that is, failed to amplify).
Identification of a sigmoidal structure is done by calculating first (FDM) and second (SDM) derivative maxima from the sigmoidal fit of the data.
In case of a sigmoidal shape, FDM is always 'right' of SDM, and in case of 'normal' amplification, the difference rarely exceeds 3 cycles
 (see 'Examples'). \code{SOD} defines 'outliers' as those with FDM < SDM or FDM - SDM > 10 or \eqn{R^2 < 0.9}.
 This function is unpublished, but works anyhow...
}

\usage{
SOD(object, remove = FALSE, verbose = TRUE, ...)
}

\arguments{
  \item{object}{an object of class 'pcrfit', 'modlist' or 'replist'.}
  \item{remove}{logical. If \code{TRUE}, the outlier runs are removed from the 'modlist'. In case of a 'replist', the data is refitted without the outlier.}
  \item{verbose}{logical. If \code{TRUE}, all analysis steps are displayed on the console.} 
  \item{...}{any other parameters to be passed to \code{\link{efficiency}} or \code{\link{replist}}.}
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

\examples{
## on a 'modlist'
## 'throw in' a failed run,
## identify and plot
dat <- reps
dat[, 10] <- rnorm(49, 1, 0.1)
ml <- modlist(dat)
res <- SOD(ml)
is.outlier(res)
plot(res, which = "single")

## remove outlier run
res <- SOD(ml, remove = TRUE)
plot(res, which = "single")

\dontrun{
## removing and updating in case
## of supplying a 'replist'
rl <- replist(ml, gl(7, 4)) 
plot(rl) # F3.1 affects the fit!
rl2 <- SOD(rl, remove = TRUE)
plot(rl2)

## to get an impression on FDM and SDM criteria,
## concatenate many different datasets and then
## calculate FDM - SDM of all runs and plot
data <- cbind(reps[1:40, ], reps2[1:40, -1], reps3[1:40, -1],
              guescini1[1:40, -1], guescini2[1:40, -1],
              rutledge[1:40, -1], batsch1[1:40, -1], 
              batsch2[1:40, -1], sisti1[1:40, -1], sisti2[1:40, -1]) 
ml <- modlist(data)
cpd1 <- sapply(ml, function(x) tryCatch(efficiency(x, plot = FALSE)$cpD1, 
                                        error = function(x) NA)) 
cpd2 <- sapply(ml, function(x) tryCatch(efficiency(x, plot = FALSE)$cpD2, 
                                        error = function(x) NA)) 

## FDM - SDM rarely exceeds 3!
barplot(cpd1 - cpd2, ylim = c(0, 10))
} 
}

\keyword{models}
\keyword{nonlinear}
