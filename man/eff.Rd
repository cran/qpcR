\name{eff}
\alias{eff}

\title{The amplification efficiency curve of a fitted object}

\description{
  Calculates the efficiency curve from the fitted model by \eqn{E = \frac{F(n)}{F(n-1)}},
  with E = efficiency, F = raw fluorescence, n = Cycle number.
}

\usage{
  eff(object, sequence = NULL, plot = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{sequence}{a sequence for the calculated curve, defaults to [min(Cycles), max(Cycles)] with 100 points per cycle.}
  \item{plot}{should the efficiency be plotted?}
 }

\value{
A list with the following components:
  \item{eff.x}{the cycle points.}
  \item{eff.y}{the efficiency curve.}
  \item{effmax}{the maximum efficiency.}
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
m <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
eff(m) #with default 100 points per cycle
eff(m, sequence = seq(5, 35, by = 0.1)) #not all data and only 10 points per cycle
}

\keyword{models}
\keyword{nonlinear}
