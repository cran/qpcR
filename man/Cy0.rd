\name{Cy0}
\alias{Cy0}

\title{Cy0 alternative to crossing points/threshold cycles as in Guescini et al}

\description{
  An alternative to the classical crossing point/threshold cycles estimation as described in Guescini et al.
  A tangent is fit to the first derivative maximum (point of inflection) of the modeled curve and its intersection 
   with the x-axis is calculated.   
}

\usage{
  Cy0(object, plot = FALSE, add = FALSE, ...)
}

\arguments{
  \item{object}{a 'drc' object.}
  \item{plot}{if \code{TRUE}, displays a plot of Cy0.}
  \item{add}{if \code{TRUE}, plot can be added to any other existing plot, i.e. as from \code{\link{pcrplot}}.}
  \item{...}{other parameters to be passed to \code{\link{pcrplot}} or \code{\link{points}}.}	
 }

\details{
  The function calculates the first derivative maximum (cpD1) of the curve and the slope and fluorescence at that point.
  Cy0 is then calculated by \eqn{Cy0 = cpD1 - (Fluo/slope)}.
}

\value{
  The Cy0 value.  
}

\author{
  Andrej-Nikolai Spiess
}

\references{
 A new real-time PCR method to overcome significant quantitative inaccuracy due to slight amplification inhibition.\cr
 Guescini M et al, \emph{BMC Bioinformatics}, 2008, \bold{9}: 326.
}
          
\examples{
### single curve with plot
m <- pcrfit(reps, 1, 2, l5())
Cy0(m, plot = TRUE)
### add to 'efficiency' plot
efficiency(m)
Cy0(m, add = TRUE)
### compare Cy0's with cpD2's
ml <- modlist(reps, fct = l5())
sapply(ml, function(x) Cy0(x))
sapply(ml, function(x) efficiency(x, plot = FALSE)$cpD2)
}

\keyword{models}
\keyword{nonlinear}
