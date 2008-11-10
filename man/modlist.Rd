\name{modlist}
\alias{modlist}

\title{Create nonlinear models from a dataframe and coerce them into a list}

\description{
 Simple function to create a list of nonlinear models from the columns of a qPCR dataframe.
 Very handy if following functions should be applied to different qPCR models, i.e. by \code{\link{sapply}}.
}

\usage{
  modlist(x, cols = NULL, fct = l4(), opt = FALSE, norm = FALSE, 
          backsub = NULL, ...)
}

\arguments{
  \item{x}{a dataframe containing the qPCR data.}
  \item{cols}{the columns (runs) to be analyzed. If \code{NULL}, all runs will be considered.}
  \item{fct}{the function used for building the model, using the function lists from the 'drc' package.} 
  \item{opt}{logical. Should model optimization take place? If \code{TRUE}, model selection is applied.} 
  \item{norm}{logical. Should the raw data be normalized to within [0, 1] before model fitting?}
  \item{backsub}{background subtraction. If \code{NULL}, not applied. Otherwise, a numeric sequence such as \code{1:10}. See 'Details'.}
  \item{...}{other parameters to be passed to \code{\link{mchoice}}.}
 }

\details{
  For a more detailed description of the functions see 'l4()' and 'pcrbatch'.
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
### subtract background using the first 8 cycles
ml <- modlist(reps, fct = l5(), backsub = 1:8)
effs <- sapply(ml, function(x) efficiency(x)$eff)
print(effs)

### 'crossing points' for the first 3 runs
### using best model from Akaike weights and normalization
ml <- modlist(reps, 2:4, fct = l4(), opt = TRUE, norm = TRUE, crit = "weights")
cps <- sapply(ml, function(x) efficiency(x)$cpD2)
print(cps)

}

\keyword{models}
\keyword{nonlinear}
