\name{modlist}
\alias{modlist}

\title{Create sigmoidal models from a dataframe and coerce them into a list}

\description{
 Simple function to create a list of nonlinear models from the columns of a qPCR dataframe.
 Very handy if following functions should be applied to different qPCR models, i.e. by \code{\link{sapply}}.
}

\usage{
  modlist(x, cyc = 1, fluo = 2:ncol(x), fct = l5())
}

\arguments{
  \item{x}{a dataframe containing the qPCR data.}
  \item{cyc}{the column containing the cycle numbers, defaults to \code{1}.}
  \item{fluo}{the column(s) containing the raw fluorescence data. If not specified, all columns will be used.}
  \item{fct}{the function used for building the model, using the function lists from the 'drc' package.}  
 }

\details{
  For a more detailed description of the functions see i.e. 'l5()'.
}

\value{
  A list with each item containing the model from each column.
}

\author{
  Andrej-Nikolai Spiess
}

\examples{
### calculate efficiencies for each run in
### the 'reps' data
ml <- modlist(reps, fct = l5())
effs <- sapply(ml, function(x) efficiency(x)$eff)
print(effs)
### 'crossing points' for the first 3 runs
### using l4 model
ml <- modlist(reps, fluo = 2:4, fct = l4())
cps <- sapply(ml, function(x) efficiency(x)$cpD2)
print(cps)
}

\keyword{IO}
\keyword{file}