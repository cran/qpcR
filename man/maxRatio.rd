\name{maxRatio}
\alias{maxRatio}

\title{The maxRatio method as in Shain et al}

\description{
The maximum ratio (MR) is determined along the interpolated curve of F(x)/F(x-1) and the corresponding
cycle number at MR is taken. A respective cycle number (FCN) is then calculated for MR.
}

\usage{
maxRatio(ml, plot = TRUE, ...)
}

\arguments{
  \item{ml}{an object of class 'modlist'.}
  \item{plot}{Should diagnostic plots be displayed?}
  \item{...}{other parameters to be passed to \code{\link{plot}}.}     
}

\details{
In the original paper the authors smooth the datapoints and then apply a cubic spline on the ratio curve, 
 in order to attain a resolution of 0.01 cycles. The function here calculates the ratio along a curve
  of the sigmoidal fit, which results in essentially the same.
}

\value{
 A list with the following components:
  \item{mr}{the maximum ratio.}
  \item{fcn}{the cycle number at \code{mr}.}
  \item{fcna}{a corrected \code{fcn}, as described in Shain et al.}
  \item{names}{the names of the runs as taken from the original dataframe.}  
}

\author{
Andrej-Nikolai Spiess
}

\references{
A new method for robust quantitative and qualitative analysis of real-time PCR.\cr
Shain & Clemens, \emph{Nucleic Acids Research}, 2008, \bold{36}, e91.
}

\examples{
ml <- modlist(reps, model = l5) 
maxRatio(ml)
}

\keyword{models}
\keyword{nonlinear}
