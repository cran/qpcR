\name{replist}
\alias{replist}

\encoding{latin1}

\title{Amalgamation of single data models into a model containing replicates}

\description{
Starting from a 'modlist' containing qPCR models from single data, \code{replist} amalgamates
 the models according to the grouping structure as defined in \code{group}. The result
 is a 'replist' with models obtained from fitting the replicates by \code{\link{pcrfit}}.
}

\usage{
replist(object, group = NULL, remove = TRUE, opt = FALSE, 
        verbose = TRUE, ...)
}

\arguments{
 \item{object}{an object of class 'modlist'.}
 \item{group}{a vector defining the replicates for each group.}
 \item{remove}{logical. If \code{TRUE}, tagged runs (that failed to fit) are automatically removed and \code{group} is updated. Should be left as is.}
 \item{opt}{should model selection be applied to the final model?}
 \item{verbose}{if \code{TRUE}, the analysis is printed to the console.}
 \item{...}{other parameters to be supplied to \code{\link{mselect}}.}    
}

\details{
As being defined by \code{group}, the 'modlist' is split into groups of runs and these amalgamated into a nonlinear model.
 If \code{remove = TRUE}, runs which had been failed to be fitted during \code{\link{modlist}} are automatically removed beforehand and \code{group} is updated
 (that is, the correpsonding entries also removed). 
 A \code{\link{nls}} model is built from the stacked raw data values. In contrast to \code{\link{curvemean}}, the fluorescence values are averaged, not the cycle numbers.
 Model selection can be applied to the final model by setting \code{opt = TRUE}.
}

\value{
An object of class 'replist' containing the replicate models of class 'nls'/'pcrfit'.  
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{modlist}}, \code{\link{pcrfit}}.
}

\examples{    
## convert 'modlist' into 'replist'
ml <- modlist(reps, model = l4)
rl <- replist(ml, group = gl(7, 4))
plot(rl)
summary(rl[[1]])

## optimize model based on Akaike weights
rl2 <- replist(ml, group = gl(7, 4), opt = TRUE, crit = "weights") 

\dontrun{
## automatic removal of unfitted runs
## Sample F3.1 is forced to fail fitting
dat <- reps
dat[, 10] <- rep(0, 49)
ml3 <- modlist(dat)
rl3 <- replist(ml3, gl(7, 4))
plot(rl3)
}
}

\keyword{models}
\keyword{nonlinear}
