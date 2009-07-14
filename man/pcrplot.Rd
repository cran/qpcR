\name{pcrplot}
\alias{pcrplot}

\title{Plotting qPCR data with fitted curves/confidence bands}

\description{
A plotting function for data of class 'pcrfit' (single curves) or 'modlist' (batch curves) displaying the
 data points, the fitted curve and (if desired) confidence/prediction bands.}

\usage{
pcrplot(object, fitted = TRUE, confband = c("none", "confidence", "prediction"),
        add = FALSE, colvec = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'pcrfit' or 'modlist'.}
  \item{fitted}{should the fitted line be displayed?}
  \item{confband}{should confidence/prediction bands be displayed? See \code{\link{confint}}.}
  \item{add}{should the curve be added to an existing plot?}
  \item{colvec}{an optional color vector for the individual curves.}   
  \item{...}{other parameters to be passed to \code{\link{plot}}, \code{\link{lines}} or \code{\link{pcrpred}}.}   
}

\value{
A qPCR plot.
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## single plot
m1 <- pcrfit(reps, 1, 2, l5)
pcrplot(m1)

## add another plot in blue
## with 99\% prediction interval
m2 <- pcrfit(reps, 1, 12, l5)
pcrplot(m2, add = TRUE, col = 4, confband = "confidence", level = 0.99)

## plot a 'modlist' batch with coloring of replicates
ml <- modlist(reps, 1, 2:13, model = l4)
pcrplot(ml, colvec = gl(3,4)) 
}

\keyword{models}
\keyword{nonlinear}
