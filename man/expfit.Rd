\name{expfit}
\alias{expfit}

\title{Calculation of PCR efficiency by fitting of an exponential model and minimization of residual variance}

\description{
  An exponential model is fit to a sliding window of defined size along the qPCR data.
  In the window with minimal residual variance (or minimal AIC), the efficiency is calculated. 
}

\usage{
  expfit(object, fitcyc = 5, plot = TRUE, crit = "resVar", start = 5, maxeff = 2, mineff = 1.5)
}

\arguments{
  \item{object}{a 'drc' object.}
  \item{fitcyc}{the size of the sliding window, default is 5.}
  \item{plot}{if TRUE the result is plotted, if FALSE the result is diplayed on the console.}	
  \item{crit}{the criterium to be minimized. Either "resVar" or "AICc".}	
  \item{start}{the cycle number to start from. Defaults to 5.} 
  \item{maxeff}{the maximum allowed efficiency of the fit.}
  \item{mineff}{the minimum allowed efficiency of the fit. See Details.}
}

\details{
The exponential growth function \eqn{f(x) = b * exp(k * x) + e} is fit to the data.
To avoid good fits in the plateau phase, the fitting procedure is stopped at cpD1 (first derivative maximum cycle).
The efficiency is calculated a) from the exponential fit with \eqn{E = exp(d)} (\code{Eff.fit}) and b) for each cycle
 within the exponential region from the raw fluorescence values by \eqn{E = \frac{F(n)}{F(n-1)}} (\code{Eff.curve}).
 The inital template fluorescence (F0) is derived from parameter b. Avoiding good fits in the noisy ground phase can be
 tweaked by increasing the \code{effmin} parameter.
}

\value{
 A list with the following components:
  \item{cyc.best}{the 'best' cycle with the minimized criterium.}
  \item{Eff.fit}{the PCR efficiency calculated from the exponential fit.}  
  \item{Eff.curve}{the PCR efficiency for each cycle within the exponential phase. See Details.}  
  \item{resVar}{the residual variance at the best cycle.}
  \item{AIC}{the AIC at the best cycle.} 
  \item{RMSE}{the root-mean-squared error at the best cycle.} 
  \item{init}{the initial template fluorescence, i.e. 'F0'.}
  \item{mod}{the exponential model from the best fit.}
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
expfit(m)
### using AIC as criterion
expfit(m, crit = "AIC")
}

\keyword{models}
\keyword{nonlinear}
