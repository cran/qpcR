\name{outlier}
\alias{outlier}

\title{Calculation of qPCR outlier cycles}

\description{
Calculates the first significant outlier cycle using the studentized residuals method.
}

\usage{
outlier(object, pval = 0.05, nsig = 3)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{pval}{the p-value for the outlier test.}
  \item{nsig}{the number of successive outlier tests. See 'Details'.}      
}

\details{
Outliers are calculated essentially as described in the reference below.
The steps are:

1) Fitting a linear model to some background cycles 1:x.\cr
2) Calculation of the studentized residuals.\cr
3) Test if the last residual is an outlier in terms of t-distribution.\cr
4) Test if the next \code{nsig} - 1 cycles are also outlier cycles.\cr
5) If so, take cycle from 3), otherwise x = x + 1  and start at 1).\cr
}

\value{
A list with the following components:
  \item{outl}{the outlier cycle.}
  \item{f.outl}{the fluorescence at \code{outl}.}    
}

\author{
Andrej-Nikolai Spiess
}

\references{
Standardized determination of real-time PCR efficiency from a single reaction set-up.
Tichopad et al., \emph{Nucleic Acids Research}, 2003, \bold{e122}.\cr       
}

\examples{
m <- pcrfit(reps, 1, 2, l5)
out <- outlier(m) 
plot(m)
abline(v = out$outl, col = 2)
abline(h = out$f.outl, col = 2)  
}

\keyword{models}
\keyword{nonlinear}
