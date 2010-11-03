\name{meltcurve}
\alias{meltcurve}

\title{Melting curve analysis and \eqn{Tm} identification}

\description{
This function conducts a melting curve analysis from the melting curve data of the real-time instrument.
The data has to be preformatted in a way that for each column of temperature values there exists a corresponding fluorescence value column.
See \code{edit(dyemelt)} for a proper format. The output is a graph displaying the raw fluorescence curve (black), the first derivative curve (red)
 and the identified melting peaks. The obtained \eqn{T_m} values are returned in a list.
}

\usage{
meltcurve(data, temps = NULL, fluos = NULL, window = NULL, 
          span = 11, ...)
}

\arguments{
  \item{data}{a dataframe containing the temperature and fluorescence data.}
  \item{temps}{a vector of column numbers reflecting the temperature values. If \code{NULL}, they are assumed to be 1, 3, 5, ... .}     
  \item{fluos}{a vector of column numbers reflecting the fluorescence values. If \code{NULL}, they are assumed to be 2, 4, 6, ... .}  	
  \item{window}{a user-defined window for the temperature region to be analyzed. See 'Details'.}
  \item{span}{the window span for peak identification. Can be tweaked to optimize \eqn{T_m} identification.}
  \item{...}{other parameters to be passed to \code{qpcR:::peaks}.}
}

\details{
The melting curve analysis in done with the following steps:\cr

1) Temperature and fluorescence values are selected in a region according to \code{window}.\cr
2) A cubic spline function (\code{\link{splinefun}}) is fit to the raw fluorescence melt values.\cr
3) The first derivative values are calculated from the spline function for each of the temperature values.\cr
4) Friedman's supersmoother (\code{\link{supsmu}}) is applied to the first derivative values.\cr
5) Melting peaks (\eqn{T_m}) values are identified by \code{qpcR:::peaks}.\cr
6) A matrix of xyy-plots is displayed using \code{qpcR:::xyy.plot}. 
}

\value{
A list with as many items as melting curves, each containing one or multiple \eqn{T_m} values.
}

\note{
The \code{peaks} function is derived from a R-Help mailing list entry in Nov 2005 by Martin Maechler.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
## default columns
data(dyemelt)
meltcurve(dyemelt, window = c(75, 87))

## selected columns
meltcurve(dyemelt, temps = c(1, 3), fluos = c(2, 4), window = c(75, 87))  
}

\keyword{models}
\keyword{nonlinear}
