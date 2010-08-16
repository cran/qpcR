\name{modlist}
\alias{modlist}

\title{Create nonlinear models from a dataframe and coerce them into a list}

\description{
Essential function to create a list of nonlinear models from the columns of a qPCR dataframe.
Very handy if following functions should be applied to different qPCR models, i.e. by \code{\link{sapply}}.
}

\usage{
modlist(x, cyc = 1, fluo = NULL, model = l4, remove = FALSE, opt = FALSE, 
        norm = FALSE, backsub = NULL, opt.method =  "LM", 
        nls.method = "port", sig.level = 0.05, crit = "ftest", ...)
}

\arguments{
  \item{x}{a dataframe containing the qPCR data or a single qPCR run of class 'pcrfit'.}
  \item{cyc}{the column containing the cycle data. Defaults to first column.}
  \item{fluo}{the column(s) (runs) to be analyzed. If \code{NULL}, all runs will be considered.}
  \item{model}{the model to be used.} 
  \item{remove}{logical. If \code{TRUE}, runs that failed to be fitted are not included in the output. See 'Details'.}
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

\details{
In case of unsuccessful model fitting and if \code{remove = FALSE} (default), the original data is included in the output, albeit with no
 fitting information. This is useful since using \code{plot.pcrfit} on the 'modlist' shows the non-fitted runs. If \code{remove = TRUE},
  the non-fitted runs are automatically removed and will thus not be displayed. See \code{\link{plot.pcrfit}} for examples.
}

\seealso{
\code{\link{pcrbatch}} for batch analysis using different methods.
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

## convert a single run to a 'modlist'
m <- pcrfit(reps, 1, 2, l5)
ml <- modlist(m)
}

\keyword{models}
\keyword{nonlinear}
