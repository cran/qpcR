\name{batchstat}
\alias{batchstat}

\encoding{latin1}

\title{Concatenating or calculating statistics on a 'pcrbatch'}

\description{
This function will either concatenate data from several \code{pcrbatch}es or calculate some user-defined statistic on the runs within a \code{pcrbatch}.
If the latter is chosen, a grouping vector must be supplied for defining the runs to be subjected to statistical analysis.
}

\usage{
batchstat(..., group = NULL,  do = c("cbind", "stat"), statfun = mean) 
}

\arguments{
  \item{...}{one or more \code{pcrbatch}es. See 'Examples'.}
  \item{group}{in case of \code{do = "stat"}, a vector defining the groups for statistical analysis.}
  \item{do}{concatenate or analyse?} 
  \item{statfun}{the statistical function to be used if \code{do = "stat"}.}   
}

\details{
\code{statfun} can be any internal R function, i.e. \code{sd}, \code{median} etc.
}

\value{
Either a concatenated dataframe (\code{do = "cbind"}), or a list containing a dataframe(s) with the statistical output for each
 factor level defined in \code{group}, if \code{do = "stat"}. 
} 

\author{
Andrej-Nikolai Spiess
}

\examples{
## create 3 'pcrbatch'es
## and concatenate
dat1 <- pcrbatch(reps, fluo = 2:5, model = l4, plot = FALSE)
dat2 <- pcrbatch(reps, fluo = 6:9, model = l4, plot = FALSE)
dat3 <- pcrbatch(reps, fluo = 10:13, model = l4, plot = FALSE)
batchstat(dat1, dat2, dat3)

## one 'pcrbatch' and doing 
## mean on replicates
## defined by 'group'
dat4 <- pcrbatch(reps, fluo = 2:9, model = l4, plot = FALSE)
GROUP <- c(1, 1, 1, 1, 2, 2, 2, 2)
batchstat(dat4, do = "stat", group = GROUP, statfun = mean)

## get the standard deviation 
batchstat(dat4, do = "stat", group = GROUP, statfun = sd)

## do stats on many 'pcrbatch'es
## All batches must have same length!
batchstat(dat1, dat2, dat3, do = "stat", 
          group = c(1, 1, 2, 2)) 
}

\keyword{models}
\keyword{nonlinear}
