\name{modlist}
\alias{modlist}

\title{Create nonlinear models from a dataframe and coerce them into a list}

\description{
Essential function to create a list of nonlinear models from the columns of a qPCR dataframe.
Very handy if following functions should be applied to different qPCR models, i.e. by \code{\link{sapply}}.
}

\usage{
modlist(x, cyc = 1, fluo = NULL, model = l4, opt = FALSE, norm = FALSE, 
        backsub = NULL, opt.method =  "LM", nls.method = "port",
        sig.level = 0.05, crit = "ftest", ...)
}

\arguments{
  \item{x}{a dataframe containing the qPCR data.}
  \item{cyc}{the column containing the cycle data. Defaults to first column.}
  \item{fluo}{the column(s) (runs) to be analyzed. If \code{NULL}, all runs will be considered.}
  \item{model}{the model to be used.} 
  \item{opt}{logical. Should model selection be applied?} 
  \item{norm}{logical. Should the raw data be normalized within [0; 1] before model fitting?}
  \item{backsub}{background subtraction. If \code{NULL}, not applied. Otherwise, a numeric sequence such as \code{1:10}. See 'Details' in \code{\link{pcrbatch}}.}
  \item{opt.method}{see \code{\link{pcrfit}}.}
  \item{nls.method}{see \code{\link{pcrfit}}.}
  \item{sig.level}{see \code{\link{mselect}}.}
  \item{crit}{see \code{\link{mselect}}.}
  \item{...}{other parameters to be passed to \code{\link{pcrfit}} or \code{\link{mselect}}.}
 }

\value{
A list with each item containing the model from each column. A 'names' item containing the column name is attached to each model.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
## calculate efficiencies for each run in
## the 'reps' data
## subtract background using the first 8 cycles
ml <- modlist(reps, model = l5, backsub = 1:8)
sapply(ml, function(x) efficiency(x, plot = FALSE)$eff)

## 'crossing points' for the first 3 runs (normalized)
##  and using best model from Akaike weights
ml <- modlist(reps, 1, 2:4, model = l5, opt = TRUE, norm = TRUE, crit = "weights" )
sapply(ml, function(x) efficiency(x, plot = FALSE)$cpD2)
}

\keyword{models}
\keyword{nonlinear}
