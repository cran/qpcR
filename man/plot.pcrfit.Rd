\name{plot.pcrfit}
\alias{plot.pcrfit}

\title{Plotting qPCR data with fitted curves/confidence bands/error bars}

\description{
A plotting function for data of class 'pcrfit' (single curves), 'modlist' (batch curves) or 'replist' (replicate curves) displaying the
 data points, the fitted curve and (if desired) confidence/prediction bands/error bars on replicates.}

\usage{
\method{plot}{pcrfit}(x, fitted = TRUE, subset = NULL, 
        confband = c("none", "confidence", "prediction"),
        errbar = c("none", "sd", "se", "conf"), add = FALSE, 
        colvec = NULL, level = 0.95, ...)
}

\arguments{
  \item{x}{an object of class 'pcrfit', 'modlist' or 'replist'.}     
  \item{fitted}{should the fitted line be displayed?}
  \item{subset}{a 2-element vector analog to \code{xlim} of the data that should be plotted.}
  \item{confband}{should confidence/prediction bands be displayed? See \code{\link{confint}}.}
  \item{errbar}{the type of error bar on the plot if replicates exist. See 'Examples'.}
  \item{add}{should the curve be added to an existing plot?}
  \item{colvec}{an optional color vector for the individual curves.} 
  \item{level}{the confidence level used for \code{confband} or \code{errbar}.}
  \item{...}{other parameters to be passed to \code{\link{plot}} or \code{\link{predict}}.}       
}

\value{
A qPCR plot. 
If \code{object} was of class 'replist', colour coding of the curves is done automatically.
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## single plot
m1 <- pcrfit(reps, 1, 2, l5)
plot(m1)

## add another plot in blue
## with 99\% prediction interval
m2 <- pcrfit(reps, 1, 12, l5)
plot(m2, add = TRUE, colvec = 4, confband = "confidence", level = 0.99)

## plot a 'modlist' batch with coloring of replicates
ml <- modlist(reps, 1, 2:13, model = l4)
plot(ml, colvec = gl(3,4)) 

## only the fitted curves
## and a subset of data
plot(ml, colvec = gl(3,4), type = "n", subset = c(10, 30))

## plot a 'replist'
ml <- modlist(reps, 1, 2:9, model = l4)
rl <- replist(ml, group = gl(2, 4))
plot(rl)

## standard deviation instead of 
## replicate points
plot(rl, type = "none", errbar = "sd")
}

\keyword{models}
\keyword{nonlinear}
