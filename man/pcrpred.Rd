\name{pcrpred}
\alias{pcrpred}

\title{Value prediction from a fitted sigmoidal model}

\description{
  After fitting the appropriate sigmoidal model, either the raw fluorescence values can be predicted from the cycle number or \emph{vice versa}. 
}

\usage{
  pcrpred(object, newdata, which = c("y", "x"), ...)

}

\arguments{
  \item{object}{a 'drc' object.}
  \item{newdata}{a vector containing the values to estimate from.}
  \item{which}{either "y" (default) for prediction of the raw fluorescence or "x" for prediction of the cycle number.}
  \item{...}{other parameters to be passed to the \code{predict.drc} function.}	
 }

\details{
Uses the \code{predict.drc} function from the 'drc' library for estimation, if the raw fluorescence values are to be estimated, and
the inversion function of either the \code{llogistic} or \code{bolztmann} models if the cycle number is to be estimated. 
}

\value{
 A vector containing the estimated values.  
}

\author{
  Andrej-Nikolai Spiess and Christian Ritz
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
### which raw fluorescence value at cycle number = 17?
pcrpred(m, 17)
### cycle numbers 20:25?
pcrpred(m, 20:25)
### which cycle at F = 4
pcrpred(m, 4, "x")
}

\keyword{models}
\keyword{nonlinear}
