\name{pcrGOF}
\alias{pcrGOF}

\title{Summarize measures for the goodness-of-fit}

\description{
  Calculates all implemented measures for the goodness-of-fit and returns them as a list.
}

\usage{
  pcrGOF(object)
}

\arguments{
  \item{object}{an object of class 'drc'.}
}

\value{
A list with the following components:
  \item{Rsq}{the R-squared value.}
  \item{AICc}{the bias-corrected Akaike Information Criterion.}
  \item{AIC}{the Akaike Information Criterion.}
  \item{resVar}{the residual variance.}
  \item{RMSE}{the root-mean-squared-error.}
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
pcrGOF(m)
}

\keyword{models}
\keyword{nonlinear}
