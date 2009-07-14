\name{sliwin}
\alias{sliwin}

\title{Calculation of PCR efficiency by the window-of-linearity method}

\description{
A linear model is fit to a sliding window of the logarithmized raw fluorescence and the regression coefficient is calculated.
At the point of maximum regression (log-linear range), the PCR efficiency is calculated. 
}

\usage{
sliwin(object, wsize = 5, border = 7, plot = TRUE)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{wsize}{the size of the sliding window, default is 5.}
  \item{border}{the -/+ border from the second derivative maximum cycle in which to do the fitting procedure.}
  \item{plot}{if TRUE the result is plotted, if FALSE the result is diplayed on the console.}	
 }

\details{
To avoid fits with a high R-squared in the baseline region, the second derivative maximum is taken as a fixpoint.
This value is consequently always near to the exponential region of the data and avoids the problem above.
The efficiency is calculated by \eqn{E = exp(slope)}, as the transformed raw data was based on the natural logarithm.
The initial template fluorescence (F0) is thus calculated by \eqn{F0 = exp(intercept)}.
}

\value{
A list with the following components:
  \item{eff}{the (maximal) PCR efficiency found within the sliding window.}
  \item{rmax}{the maximum R-squared.}
  \item{init}{the initial template fluorescence F0.}  
}

\author{
Andrej-Nikolai Spiess
}

\references{
Assumption-free analysis of quantitative real-time polymerase chain reaction (PCR) data.
Ramakers C et al., \emph{Neuroscience Letters}, 2003, \bold{339}, 62-66.
}

\examples{
m <- pcrfit(reps, 1, 2, l5)
sliwin(m)
}

\keyword{models}
\keyword{nonlinear}
