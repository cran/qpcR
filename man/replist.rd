\name{replist}
\alias{replist}

\encoding{latin1}

\title{Amalgamation of single data models to a model containing replicates}

\description{
Starting from a 'modlist' containing qPCR models from single data, \code{replist} amalgamates
 the models according to the grouping structure as defined in \code{group}. The result
 is a 'modlist' with models obtained from fitting the replicates by \code{\link{pcrfit}}.
}

\usage{
replist(object, group)
}

\arguments{
\item{object}{an object of class 'modlist'.}
\item{group}{a vector defining the replicates for each group.}
}

\details{
As being defined by \code{group}, the raw data of the curves are averaged and starting values for the averaged values are
 calculated by \code{\link{optim}}. Finally, a \code{\link{nls}} model is built from the \code{\link{stack}}ed raw data values using
 these starting values.
}

\value{
An object of class 'modlist' containing the replicate models of class 'nls'/'pcrfit'.  
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{modlist}}, \code{\link{pcrfit}}.
}

\examples{
\dontrun{
ml <- modlist(reps, model = l4)
rl <- replist(ml, group = gl(7, 4))
plot(rl)
summary(rl[[1]])
}
}

\keyword{models}
\keyword{nonlinear}
