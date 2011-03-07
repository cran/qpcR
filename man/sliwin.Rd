\name{sliwin}
\alias{sliwin}

\title{Calculation of qPCR efficiency by the 'window-of-linearity' method}

\description{
A linear model is fit to a sliding window of defined size(s) and within a defined border around a fixpoint cycle. Regression coefficients are calculated for each window, and at the point of maximum regression (log-linear range) or least variation in slope, parameters such as PCR efficiency and initial template fluorescence are calculated. From version 1.3-5, an approach "not unlike" to Ruijter et al. (2009) has been implemented, in which baseline values can be iteratively subtracted from the data prior to fitting the sliding window. See 'Details' for more information.  
}

\usage{
sliwin(object, wsize = 6, fix = c("top", "cpD2"), type = c("rsq", "slope"),
       border = c(0, 8), base = 0, plot = TRUE, ...)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{wsize}{the size(s) of the sliding window(s), default is \code{6}. A sequence such as \code{4:6} can be used to optimize the window size.}
  \item{fix}{which method to use for the border fixpoint. If \code{top}, it is the take-off point, if \code{cpD2}, the second derivative maximum. See 'Details'.}
  \item{type}{selection of the window with best baseline + maximum R-square (\code{rsq}) or best baseline + minimal variance in slope + maximum R-square (\code{"slope"}).}
  \item{border}{the -/+ cycle border from the \code{fix}point in which to do the fitting procedure. Can have significant impact on the efficiency outcome, see 'Details'.}
  \item{base}{either \code{0} for no baseline optimization, or a scalar defining multiples of the standard deviation of all baseline points \code{1:top (take-off point)}, which are iteratively subtracted from the raw data. See 'Details' and 'Examples'.}
  \item{plot}{if \code{TRUE}, the result is plotted with the logarithmized curve, sliding window, regression line and baseline.}	
  \item{...}{only used internally for passing the parameter matrix.}
}

\details{
To avoid fits with a high \eqn{R^2} in the baseline region, some fixpoint in the data must be taken. In case of \code{fix = "top"} (default), this is the take-off point from function \code{\link{takeoff}}. This value is always at the start of the exponential region, while the second derivative maximum with \code{fix = "cpD2"} is usually at the end. In the latter case, \code{border} values such as \code{c(-8, 0)} make sense. The efficiency is calculated by \eqn{E = 10^{slope}}, as the transformed raw data are based on \code{log10}. The initial template fluorescence (F0) is thus calculated by \eqn{F0 = 10^{intercept}}. For the baseline optimization, 100 baseline values \eqn{Fb_i} are interpolated in the range of the data:
\deqn{F_{min} \le Fb_i \le base \cdot \sigma(F_1...F_{top})}
For each iteration, the best window in terms of linear regression is found using the window sizes in \code{wsize} and the the borders in \code{border}, if \code{type = "rsq"}. If \code{type = "slope"}, the baseline value delivering the smallest variance in the slope of all sliding windows is selected, and then the window with highest \eqn{R^2}. This approach is quite similar to the one in Ruijter et al. (2009) but has to be tweaked in order to obtain the same values as in the 'LinRegPCR' software. Especially the \code{border} value has significant influence on the calculation of the best window's efficiency value. 
}

\value{
A list with the following components:
  \item{eff}{the optimized PCR efficiency found within the sliding window.}
  \item{rsq}{the maximum R-squared.}
  \item{init}{the initial template fluorescence F0.}  
  \item{base}{the optimized baseline value.}
  \item{window}{the best window found within the \code{border}s.}  
  \item{parMat}{a matrix containing the parameters as above for each iteration.}
}

\author{
Andrej-Nikolai Spiess
}

\references{
Assumption-free analysis of quantitative real-time polymerase chain reaction (PCR) data.\cr
Ramakers C, Ruijter JM, Deprez RH, Moorman AF.\cr
\emph{Neurosci Lett} (2003), \bold{339}: 62-65.\cr

Amplification efficiency: linking baseline and bias in the analysis of quantitative PCR data.\cr
Ruijter JM, Ramakers C, Hoogaars WM, Karlen Y, Bakker O, van den Hoff MJ, Moorman AF.\cr
\emph{Nucleic Acids Res} (2009),  \bold{37}: e45
}

\examples{
## sliding window of size 5
## between take-off point and 10 cycles 
## upstream, no baseline optimization
m1 <- pcrfit(reps, 1, 2, l4)
sliwin(m1, wsize = 5, border = c(0, 10), base = 0)

\dontrun{
## optimizing with window sizes of 4 to 6,
## between 0/+8 from take-off point,
## and baseline up to 2 standard deviations
sliwin(m1, wsize = 4:6, border = c(0, 8), 
       base = 2)
}
}

\keyword{models}
\keyword{nonlinear}
