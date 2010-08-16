\name{replist}
\alias{replist}

\encoding{latin1}

\title{Amalgamation of single data models to a model containing replicates}

\description{
Starting from a 'modlist' containing qPCR models from single data, \code{replist} amalgamates
 the models according to the grouping structure as defined in \code{group}. The result
 is a 'replist' with models obtained from fitting the replicates by \code{\link{pcrfit}}.
}

\usage{
replist(object, group = NULL, opt = FALSE, verbose = TRUE, ...)
}

\arguments{
\item{object}{an object of class 'modlist'.}
\item{group}{a vector defining the replicates for each group.}
\item{opt}{should model selection be applied to the final model?}
\item{verbose}{if \code{TRUE}, the analysis is printed to the console.}
\item{...}{other parameters to be supplied to \code{\link{mselect}}.}

}

\details{
As being defined by \code{group}, the raw data of the curves are averaged and starting values for the averaged values are
 calculated by \code{\link{optim}}. Finally, a \code{\link{nls}} model is built from the \code{\link{stack}}ed raw data values using
 these starting values. In contrast to \code{\link{curvemean}}, the fluorescence values are averaged, not the cycle numbers.
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
ml <- modlist(reps, model = l4)
rl <- replist(ml, group = gl(7, 4))
plot(rl)
summary(rl[[1]])

## optimize model based on Akaike weights
rl2 <- replist(ml, group = gl(7, 4), opt = TRUE, crit = "weights")

}

\keyword{models}
\keyword{nonlinear}
