\name{pcropt2}
\alias{pcropt2}

\title{Random elimination and sigmoidal re-fitting of qPCR data (qPCR jackknifing)}

\description{
From a fitted model, \code{nsample} number of random datapoints are eliminated \code{B} number of times, 
 and the sigmoidal model is refit on the reduced dataset. The usual features (efficiency, crossing points, goodness of fit etc.)
 are calculated on the models and summarized at the end by mean/standard deviaton.
A plot containing boxplots for all features including a confidence interval (red lines) can be displayed. 
}

\usage{
pcropt2(object, nsample = 5, B = 100, plot = TRUE, alpha = 0.05, ...)
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{nsample}{numeric. The number of datapoints to remove in each iteration.} 
  \item{B}{numeric. The number of jackknifes.} 
  \item{plot}{logical. If \code{TRUE}, boxplots of the resulting values are displayed.}
  \item{alpha}{the alpha for the confidence interval, defaults to 95\%.}
  \item{...}{other parameters to be passed on to \code{\link{efficiency}}, \code{\link{pcrfit}} or \code{\link{boxplot}}.}
}

\details{
It has been shown by Rutledge (2004) that the estimation of PCR efficiency gives more
 realistic values when the number of plateau cycles are decreased. The function here goes even further by eliminating cycles at all
 kinds of positions and refitting the sigmoidal model. From a theoretical point of view, the jackknifing should give more realistic estimates
 of the evaluated features and their variance. 
}

\value{
A list containing the following items:
  \item{raw}{the calulated features in each iteration}.
  \item{mean}{the mean values from \code{raw}.}
  \item{sd}{the standard deviation from \code{raw}.}
  \item{upper}{the upper confidence interval for \code{raw}.}  
  \item{lower}{the lower confidence interval for \code{raw}.}  
}

\author{
Andrej-Nikolai Spiess
}

\references{
Jackknife, bootstrap and other resampling methods in regression analysis.\cr
C.F.J. Wu\cr
\emph{Ann. Stat}. (1986), \bold{14}: 1261-1295.
}



\examples{
## Removing 5 datapoints and 
## 200 iterations
m <- pcrfit(reps, 1, 2, l4)
pcropt2(m, nsample = 5, B = 200)
}


\keyword{models}
\keyword{nonlinear}
