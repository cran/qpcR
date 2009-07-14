\name{expcomp}
\alias{expcomp}

\title{Comparison of all available sigmodal models by RMSE within the exponential region}

\description{
The exponential region of the qPCR data is identified by the studentized outlier method, as in \link{expfit}.
The root-mean-squared-error of all available sigmoidal models within this region is then calculated.
The result of the fits are plotted and models returned in order of ascending RMSE.
}

\usage{
expcomp(object, ...)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{...}{other parameters to be passed to \code{expfit}.} 
}

\details{
The following sigmoidal models are fitted:
b3, b4, b5, l3, l4, l5, w3, w4, baro5
}

\value{
A dataframe with names of the models, in ascending order of RMSE.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
m <- pcrfit(reps, 1, 2, l4)
expcomp(m)
}

\keyword{models}
\keyword{nonlinear}
