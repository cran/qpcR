\name{modlist}
\alias{modlist}

\title{Create nonlinear models from a dataframe and coerce them into a list}

\description{
 Simple function to create a list of nonlinear models from the columns of a qPCR dataframe.
 Very handy if following functions should be applied to different qPCR models, i.e. by \code{\link{sapply}}.
}

\usage{
  modlist(x, cyc = 1, fluo = 2:ncol(x), fct = l5(), opt = FALSE, norm = FALSE)
}

\arguments{
  \item{x}{a dataframe containing the qPCR data.}
  \item{cyc}{the column containing the cycle numbers, defaults to \code{1}.}
  \item{fluo}{the column(s) containing the raw fluorescence data. If not specified, all columns will be used.}
  \item{fct}{the function used for building the model, using the function lists from the 'drc' package.} 
  \item{opt}{logical. Should model optimization take place? If \code{TRUE}, model selection is applied.} 
  \item{norm}{logical. Should the raw data be normalized to within [0, 1] before model fitting?}
 }

\details{
  For a more detailed description of the functions see i.e. 'l5()'.
}

\value{
  A list with each item containing the model from each column. A 'names' attribute containing the column name is attached to each model.
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
### using best model and normalization
ml <- modlist(reps, fluo = 2:4, fct = l4(), opt = TRUE, norm = TRUE)
cps <- sapply(ml, function(x) efficiency(x)$cpD2)
print(cps)

}

\keyword{IO}
\keyword{file}
