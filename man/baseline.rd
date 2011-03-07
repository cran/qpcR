\name{baseline}
\alias{baseline}

\title{Baselining and refitting six-parameter sigmoidal models}

\description{
qPCR run(s) that were fit with either the \code{\link{l6}} or \code{\link{b6}} six-parameter models offer the possibility of baselining,
 that is, subtracting the offset (parameter \code{c}) and slope (parameter \code{k}) from the response (fluorescence) values.
 This is supposed to model the baseline region of the PCR curves better than other models in which a baseline slope is missing.
 The function works on individual models of class 'pcrfit' or on 'modlist's. If \code{refit = TRUE} (default), the baselined data
 is used to refit the model.   
}

\usage{
baseline(object, refit = TRUE, refit.model = NULL, verbose = TRUE, ...)
}

\arguments{
  \item{object}{an object of class 'pcrfit' or 'modlist'.}
  \item{refit}{logical. If \code{TRUE}, the baselined data is used for refitting the model.}
  \item{refit.model}{an optional new model used for refitting. If \code{NULL}, the original model, i.e. \code{l6} or \code{b6} from the object is used.}
  \item{verbose}{logical. If \code{TRUE}, the analysis steps are displayed in the console window.} 
  \item{...}{other parameters to be passed to \code{\link{pcrfit}}.} 
}

\details{
This function is experimental and ongoing work will show if it is feasible.
}

\value{
Either a (refitted) model or a 'modlist' containing the (refitted) models.
The returned object or 'modlist' items have a new \code{$DATA.base} object containing the baselined data. See 'Examples'.
}

\author{
Andrej-Nikolai Spiess, based on an idea from Eric Shain.
}

\examples{
## baselining a single model
m <- pcrfit(reps, 1, 2, l6)
res <- baseline(m)
res$DATA.base

## baselining a 'modlist'
ml <- modlist(reps, fluo = 2:5, model = l6)
res <- baseline(ml)
plot(res) 
}

\keyword{models}
\keyword{nonlinear}
