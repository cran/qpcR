\name{expfit}
\alias{expfit}

\title{Calculation of PCR efficiency by fitting an exponential model}

\description{
  An exponential model is fit to a window of defined size on the qPCR raw data.
  The window is identified either by the 'studentized outlier' method as described in Tichopad et al. (2003),
 the 'midpoint' method (Peirson et al., 2003) or by subtracting the difference of cpD1 and cpD2 from cpD2 ('ERBCP', unpublished).
}

\usage{
expfit(object, method = c("outlier", "midpoint", "ERBCP"), pval = 0.05, 
n.outl = 3, n.ground = 10, corfact = 1, fix = c("top", "bottom", "middle"), 
nfit = 5, plot = TRUE, ...) 
}

\arguments{
  \item{object}{a 'drc' object.}
  \item{method}{one of the three possible methods to be used for defining the position of the fitting window.}
  \item{pval}{for method "outlier", the p-value for the outlier test.}
  \item{n.outl}{for method "outlier", the number of successive outlier cycles.}
  \item{n.ground}{for method "midpoint", the number of cycles in the noisy ground phase to calculate the standard deviation from.}
  \item{corfact}{for method "ERBCP", the correction factor for finding the exponential region. See 'Details'.}
  \item{fix}{for methods 'midpoint' and 'ERBCP', the orientation of the fitting window based on the identified point. See 'Details'.} 
  \item{nfit}{the size of the fitting window.}
  \item{plot}{logical. If \code{TRUE}, a graphical display of the curve and the fitted region is shown.}
  \item{...}{other parameters to be passed to the plotting function.} 
}

\details{
The exponential growth function \eqn{f(x) = b * exp(k * x) + e} is fit to the data.
Calls \code{outlier} and \code{studres} for calculation of the studentized residuals and 'outlier' cycles, and \code{midpoint}
 for calculation of the exponential phase 'midpoint'. For method 'ERBCP' (Exponential Region By Crossing Points), the exponential region
 is calculated by \eqn{expR = cpD2 - \code{corfact} * (cpD1-cpD2)}. The efficiency is calculated a) from the exponential fit with \eqn{Eff = exp(k)} and b) for each 
 cycle within the exponential region from the raw fluorescence values by \eqn{Eff = \frac{F(n)}{F(n-1)}}. 
 The inital template fluorescence (F0) is derived from parameter b.
}

\value{
 A list with the following components:
  \item{point}{the point within the exponential region as identified by one of the three methods.}
  \item{cycles}{the cycles of the identified region as defined by \code{method}, \code{fix} and \code{nfit}.}
  \item{eff}{the efficiency calculated from the exponential fit.}
  \item{eff.cycles}{the efficiency of all points within the identified region.}
  \item{AIC}{the Akaike Information Criterion of the fit.}
  \item{resVar}{the residual variance of the fit.}
  \item{RMSE}{the root-mean-squared-error of the fit.}
  \item{init}{the initial template fluorescence.}
  \item{mod}{the exponential model.}
}

\author{
  Andrej-Nikolai Spiess
}

\references{
  Standardized determination of real-time PCR efficiency from a single reaction set-up.
  Tichopad et al., \emph{Nucleic Acids Research}, 2003, \bold{e122}.\cr
  Experimental validation of novel and conventional approaches to quantitative real-time PCR data analysis.
  Peirson et al., \emph{Nucleic Acids Research}, 2003, \bold{e73}.
}

\examples{
m <- pcrfit(reps, 1, 2, l5())
### using 'outlier' method
expfit(m)
### 'midpoint' method and 7 cycle window
expfit(m, method = "midpoint", nfit = 7)
### 'ERBCP' method with window centered around
### fixpoint
expfit(m, method = "ERBCP", fix = "middle")
}

\keyword{models}
\keyword{nonlinear}
