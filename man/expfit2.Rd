\name{expfit2}
\alias{expfit2}

\title{Calculation of PCR efficiency by fitting of an exponential model after outlier cycle detection}

\description{
  The outlier cycles (and thus the region when the PCR curve rises significantly above background) are identified by 
  the studentized residual method as in Tichopad et al. (2003).
  An exponential model is then fit to the next defined number of cycles. 
}

\usage{
  expfit2(object, fitcyc = 5, plot = TRUE, pval = 0.01, ...)
}

\arguments{
  \item{object}{a 'drc' object.}
  \item{fitcyc}{the number of cycles used for exponential fitting.}
  \item{plot}{if TRUE the result is plotted, if FALSE the result is diplayed on the console.}
  \item{pval}{the p-value for the outlier significance of the studentized residuals.}	
  \item{...}{other arguments to be passed to \code{multdrc}.}  
 }

\details{
The exponential growth function f(x) = b * exp(d * x) + e is fit to the data.
Calls \code{outlier} and \code{studres} for calculation of the studentized residuals and definition of the outlier cycles.
The efficiency is calculated a) from the exponential fit with \eqn{E = exp(d)} (\code{Eff.fit}) and b) for each cycle within
 the exponential region from the raw fluorescence values by \eqn{E = \frac{F(n)}{F(n-1)}} (\code{Eff.curve}). The inital
 template fluorescence (F0) is derived from parameter b. If the detection of outlier cycles is problematic (i.e. due to high noise), one can improve results by lowering the p-value.
}

\value{
 A list with the following components:
  \item{outlier}{the outlier cycle.}
  \item{Eff.fit}{the PCR efficiency calculated from the exponential fit.}  
  \item{Eff.curve}{the PCR efficiency for each cycle within the exponential phase. See Details.}  
  \item{resVar}{the residual variance of the fit.}  
  \item{AIC}{the Akaike Information Criterion of the fit.} 
  \item{RMSE}{the root-mean-squared-error of the fit.} 
  \item{init}{the initial template fluorescence, i.e. 'F0'.}
  \item{mod}{the exponential model from the best fit.}
}

\author{
  Andrej-Nikolai Spiess
}

\references{
  Standardized determination of real-time PCR efficiency from a single reaction set-up.
  Tichopad et al., \emph{Nucleic Acids Research}, 2003, \bold{e122}.
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
expfit2(m)
}

\keyword{models}
\keyword{nonlinear}
