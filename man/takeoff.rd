\name{takeoff}
\alias{takeoff}

\title{Calculation of the qPCR takeoff point}

\description{
Calculates the first significant cycle of the exponential region (takeoff point) using the studentized residuals method from Tichopad et al (2003).
}

\usage{
takeoff(object, pval = 0.05, nsig = 3)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{pval}{the p-value for the takeoff test.}
  \item{nsig}{the number of successive takeoff tests. See 'Details'.}      
}

\details{
Takeoff points are calculated essentially as described in the reference below.
The steps are:

1) Fitting a linear model to some background cycles 1:x.\cr
2) Calculation of the studentized residuals.\cr
3) Test if the last residual is an outlier in terms of t-distribution.\cr
4) Test if the next \code{nsig} - 1 cycles are also outlier cycles.\cr
5) If so, take cycle from 3), otherwise x = x + 1  and start at 1).\cr
}

\value{
A list with the following components:
  \item{top}{the takeoff point.}
  \item{f.top}{the fluorescence at \code{top}.}    
}

\author{
Andrej-Nikolai Spiess
}

\references{
Standardized determination of real-time PCR efficiency from a single reaction set-up.
Tichopad et al., \emph{Nucleic Acids Research}, 2003, \bold{e122}.\cr       
}

\examples{
m1 <- pcrfit(reps, 1, 2, l5)
res1 <- takeoff(m1) 
plot(m1)
abline(v = res1$top, col = 2)
abline(h = res1$f.top, col = 2)  
}

\keyword{models}
\keyword{nonlinear}
