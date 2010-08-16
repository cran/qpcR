\name{rep2mod}
\alias{rep2mod}

\encoding{latin1}

\title{Converts a 'replist' back to a 'modlist'}

\description{
This function is essentially the reverse of \code{\link{modlist}}. An object of class 'replist' is converted back to a 'modlist',
 with the original grouping structure attached. 
}

\usage{
rep2mod(rl)
}

\arguments{
\item{rl}{a 'replist' containing replicates.}
}

\details{
The returned 'modlist' has an attribute 'group' attached, which is a vector containing the original grouping 
 (number of groups, number of replicates). This can be queried by \code{attr(name, "group")}.
}

\value{
A 'modlist' containing all single runs.
}

\author{
Andrej-Nikolai Spiess
}


\examples{
## create a 'modlist'
ml <- modlist(reps, 1, 2:5, l4)

## convert into 'replist'
rl <- replist(ml, group = rep(1, 4))
par(mfrow = c(2, 1))
plot(rl, main = "A 'replist'")

## convert back to a 'modlist'
ml2 <- rep2mod(rl) 
plot(ml2, main = "A 'modlist'")
}

\keyword{models}
\keyword{nonlinear}
