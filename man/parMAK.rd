\name{parMAK}
\alias{parMAK}

\title{Parameters that can be changed to tweak the mak2/mak3 methods}

\description{
A function to set (at the moment) three parameters that change the selfStart performance of mak2/mak3 models. This is feasible if the fitting fails to converge, or if there are other problems such as a high background noise level.
}

\usage{
parMAK(SS.offset = 0, SS.method = "LM", SS.deriv = c("sigfit", "spline")) 
}

\arguments{
  \item{SS.offset}{number of cycles to add to the second derivative maximum cut-off point.}
  \item{SS.method}{fitting algorithm to be used for selfStart parameter estimation. Either \code{"LM"} or any of the methods in \code{\link{optim}}.}  
  \item{SS.deriv}{method to calculate the second derivative maximum cut-off point. Either \code{"sigfit"} based on 4-par sigmoidal fit or \code{"spline"} based on smoothing and estimation from the smoothed data by two-times \code{\link{diff}}.}  
}

\details{
Calling the function writes a list named \code{parMAKs} with the three parameters to the global environment. In the \bold{mak2} and \bold{mak3} methods, these parameters are passed to \code{$ssFct} and \code{$fct_ssFct}. Selecting \code{SS.deriv = "spline"} usually leads to 1-2 cycles less for cutoff compared to \code{"sigfit"}. The user should try which is the best setting by using some calibration data and choosing the method delivering the highest \eqn{R^2} in the linear regression plot.
}

\value{
A list with three items is written to the global environment.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
## fitting a mak2 model (no slope parameter)
m1 <- pcrfit(reps, 1, 2, mak2)
plot(m1)
BIC(m1)

## fitting a mak3 model (with slope parameter)
m2 <- pcrfit(reps, 1, 2, mak3)
plot(m2, add = TRUE, col = 2)
BIC(m2)

## using a spline fit
parMAK(SS.deriv = "spline")
m3 <- pcrfit(reps, 1, 2, mak3)
}

\keyword{models}
\keyword{nonlinear}
