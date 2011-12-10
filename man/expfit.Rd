\name{expfit}
\alias{expfit}

\title{Calculation of PCR efficiency by fitting an exponential model}

\description{
An exponential model is fit to a window of defined size on the qPCR raw data.
The window is identified either by the 'studentized outlier' method as described in Tichopad \emph{et al}. (2003), the 'midpoint' method (Peirson \emph{et al}., 2003) or by subtracting the difference of cpD1 and cpD2 from cpD2 ('ERBCP', unpublished).
}

\usage{
expfit(object, method = c("outlier", "midpoint", "ERBCP"), pval = 0.05, 
       n.outl = 3, n.ground = 1:5, corfact = 1, 
       fix = c("top", "bottom", "middle"), nfit = 5, plot = TRUE, ...) 
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{method}{one of the three possible methods to be used for defining the position of the fitting window.}
  \item{pval}{for \code{method = "outlier"}, the p-value for the outlier test.}
  \item{n.outl}{for \code{method = "outlier"}, the number of successive outlier cycles.}
  \item{n.ground}{for \code{method = "midpoint"}, the number of cycles in the noisy ground phase to calculate the standard deviation from.}
  \item{corfact}{for \code{method = "ERBCP"}, the correction factor for finding the exponential region. See 'Details'.}
  \item{fix}{for methods "midpoint" and "ERBCP", the orientation of the fitting window based on the identified point. See 'Details'.} 
  \item{nfit}{the size of the fitting window.}
  \item{plot}{logical. If \code{TRUE}, a graphical display of the curve and the fitted region is shown.}
  \item{...}{other parameters to be passed to the plotting function.} 
}

\details{
The exponential growth function \eqn{f(x) = b \cdot exp(k \cdot x) + e} is fit to the data. Calls \code{\link{takeoff}} for the calculation of the studentized residuals and 'outlier' cycle, and \code{midpoint} for calculation of the exponential phase 'midpoint'. For method 'ERBCP' (Exponential Region By Crossing Points), the exponential region is calculated by \eqn{expR = cpD2 - \code{corfact} \cdot (cpD1-cpD2)}. The efficiency is calculated a) from the exponential fit with \eqn{E = exp(k)} and b) for each cycle within the exponential region from the raw fluorescence values by \eqn{E_n = \frac{F_n}{F_{n-1}}}. The inital template fluorescence \eqn{F_0} is derived from \eqn{F_0 = b \cdot exp(k)}.
}

\value{
 A list with the following components:
  \item{point}{the point within the exponential region as identified by one of the three methods.}
  \item{cycles}{the cycles of the identified region as defined by \code{method}, \code{fix} and \code{nfit}.}
  \item{eff}{the efficiency calculated from the exponential fit.}
  \item{eff.cycles}{the efficiencies of all points within the identified region.}
  \item{AIC}{the Akaike Information Criterion of the fit.}
  \item{resVar}{the residual variance of the fit.}
  \item{RMSE}{the root-mean-squared-error of the fit.}
  \item{init}{the initial template fluorescence.}
  \item{mod}{the exponential model of class 'nls'.}
}

\author{
Andrej-Nikolai Spiess
}

\references{
Standardized determination of real-time PCR efficiency from a single reaction set-up.\cr
Tichopad A, Dilger M, Schwarz G & Pfaffl MW.\cr
\emph{Nucleic Acids Research} (2003), \bold{31}:e122.

Comprehensive algorithm for quantitative real-time polymerase chain reaction.\cr
Zhao S & Fernald RD.\cr
\emph{J Comput Biol} (2005), \bold{12}:1047-64.
}

\examples{
## using 'outlier' method
m1 <- pcrfit(reps, 1, 2, l5)
expfit(m1)

## 'midpoint' method and 7 cycle window
expfit(m1, method = "midpoint", nfit = 7)

## 'ERBCP' method with window centered around
## fixpoint
expfit(m1, method = "ERBCP", fix = "middle")
}

\keyword{models}
\keyword{nonlinear}
