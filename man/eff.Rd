\name{eff}
\alias{eff}

\title{The amplification efficiency curve of a fitted object}

\description{
Calculates the efficiency curve from the fitted object by \eqn{E = \frac{F(n)}{F(n-1)}},
 with E = efficiency, F = raw fluorescence, n = Cycle number.
}

\usage{
eff(object, sequence = NULL, plot = FALSE)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{sequence}{a 3-element vector (from; to; by) defining the sequence for the efficiency curve, defaults to [min(Cycles), max(Cycles)] with 100 points per cycle.}
  \item{plot}{should the efficiency be plotted?}
 }

\value{
A list with the following components:
  \item{eff.x}{the cycle points.}
  \item{eff.y}{the efficiency values at \code{eff.x}.}
  \item{effmax.x}{the cycle number with the highest efficiency.}
  \item{effmax.y}{the maximum efficiency.}   
}

\author{
Andrej-Nikolai Spiess
}      

\examples{
m <- pcrfit(reps, 1, 2, l5)

## with default 100 points per cycle
eff(m, plot = TRUE) 

## not all data and only 10 points per cycle
eff(m, sequence = c(5, 35, 0.1), plot = TRUE) 
}

\keyword{models}
\keyword{nonlinear}
