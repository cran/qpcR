\name{midpoint}
\alias{midpoint}

\title{Calculation of the 'midpoint' region according to Peirson et al. (2003)}

\description{
Calculates the exponential region midpoint using the algorithm described under 'References'.
}

\usage{
midpoint(object, noise.cyc = 1:5)
}

\arguments{
  \item{object}{a fitted object.}
  \item{noise.cyc}{the cycles defining the background noise.}      
}

\details{
The 'midpoint' region is calculated by \deqn{F_{noise} \times \sqrt{\frac{F_{max}}{F_{noise}}}}
 with Fnoise = the standard deviation of the background cycles and Fmax = the maximal fluorescence.
}

\value{
 A list with the following components:
  \item{f.mp}{the 'midpoint' fluorescence.}
  \item{cyc.mp}{the 'midpoint' cycle, as predicted from \code{f.mp}.}    
}

\author{
Andrej-Nikolai Spiess
}

\references{
Experimental validation of novel and conventional approaches to quantitative real-time PCR data analysis.
Peirson et al., \emph{Nucleic Acids Research}, 2003, \bold{e73}.  
}

\examples{
m1 <- pcrfit(reps, 1, 2, l5)
mp <- midpoint(m1) 
plot(m1)
abline(h = mp$f.mp, col = 2)
abline(v = mp$mp, col = 2)  
}

\keyword{models}
\keyword{nonlinear}
