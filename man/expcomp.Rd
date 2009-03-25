\name{expcomp}
\alias{expcomp}

\title{Comparison of several sigmodal models by RMSE within the exponential region}

\description{
  The exponential region of the qPCR data is identified by the studentized outlier method, as in \link{expfit}.
  The root-mean-squared-errors of 14 different sigmoidal models within this region are then calculated.
  The result of the fits are plotted and the models are returned in order of ascending RMSE.
}

\usage{
  expcomp(object, ...)
}

\arguments{
  \item{object}{a 'drc' object.} 
  \item{...}{other parameters to be passed to \code{expfit}.} 
}

\details{
The following sigmoidal models are fitted:
b3(), b4(), b5(), l3(), l4(), l5(), LL2.3(), LL2.4(), LL2.5(), baro5(), w3(), w4(), W2.3(), W2.4()
}

\value{
A dataframe with names of the models and the respective RMSEs.
}

\author{
  Andrej-Nikolai Spiess
}

\examples{
m <- pcrfit(reps, 1, 2, l5())
expcomp(m)
}

\keyword{models}
\keyword{nonlinear}
