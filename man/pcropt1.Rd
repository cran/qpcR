\name{pcropt1}
\alias{pcropt1}

\title{Combinatorial elimination of plateau and ground phase cycles}

\description{
  The estimation of PCR efficiency and calculation of initial fluorescence (F0) is analyzed by refitting the (optimized) model on subsets of the data,
 thereby using all possible combinations of datapoints. The estimated parameters are then collated in a dataframe.
 This is intended to be the prerequisite for finding the optimal datapoints that minimize the fit or exhibit the best correlation to a calibration curve.
}

\usage{
  pcropt1(object, fact = 3, opt = FALSE, ...)
}

\arguments{
  \item{object}{a 'drc' object.}
  \item{fact}{numeric. The multiplier for the scan border. See 'Details'.} 
  \item{opt}{logical. If \code{true}, model selection is applied for each combination of cycles. Beware: Slow!!} 
  \item{...}{other parameters to be passed on to \code{\link{efficiency}}, i.e. where to take the PCR efficiency from.}
}

\details{
  It has been shown by Rutledge (2004) that the estimation of PCR efficiency gives more
  realistic values when the number of plateau cycles are decreased. This paradigm is the basis for
  this function, but we also consider the cycles in the ground phase and all combinations between ground/plateau cycles.
  All datapoints between the lower border cpD1 - \code{fact} * (cpD1 - cpD2) and upper border cpD1 + \code{fact} * (cpD1 - cpD2) are cycled through.
}

\value{
  A dataframe with the border values, the AICc, the residual variance, the efficiency and the estimated F(0) from the exponential and sigmoidal model.
}

\author{
Andrej-Nikolai Spiess
}

\references{
  Sigmoidal curve fitting redefines quantitative real-time PCR with the prospective of developing automated high-throughput applications. 
  Rutledge RG, \emph{Nucleic Acids Research}, 2004, \bold{e178}.
}

\seealso{
  The function \code{\link{efficiency}} that is called by this function.
}

\examples{
\dontrun{
### Using one model throughout
m <- pcrfit(reps, 1, 2, l5())
temp <- pcropt1(m)
edit(temp)
### Selecting the best model for each combination
m <- pcrfit(reps, 1, 2, l4())
temp <- pcropt1(m, opt = TRUE)
edit(temp)
}
}

\keyword{models}
\keyword{nonlinear}
