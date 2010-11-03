\name{meanlist}
\alias{meanlist}

\encoding{latin1}

\title{Amalgamation of single data models into an averaged model}

\description{
Starting from a 'modlist' containing qPCR models from single data, \code{meanlist} amalgamates
 the models according to the grouping structure as defined in \code{group}. The result
 is a 'modlist' with models obtained from averaging the replicates by \code{\link{pcrfit}}.
}

\usage{
meanlist(object, group, type = c("mean", "median"))
}

\arguments{
\item{object}{an object of class 'modlist'.}
\item{group}{a vector defining the replicates for each group.}
\item{type}{how to average the data.}
}

\details{
As being defined by \code{group}, the average data of the curves is subjected to \code{\link{pcrfit}} and a new modlist with the
 averaged models is created. Similar to \code{\link{replist}} but does not contain the replicates within the 'nls' model but the averaged
  model with only ONE curve.
}

\value{
An object of class 'modlist' containing the averaged models of class 'nls'/'pcrfit'.  
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{modlist}}, \code{\link{replist}}.
}

\examples{
ml <- modlist(reps, model = l4)
res <- meanlist(ml, group = gl(7, 4))
plot(res)
efficiency(res[[1]])
}

\keyword{models}
\keyword{nonlinear}
