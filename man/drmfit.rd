\name{drmfit}
\alias{drmfit}

\title{The drm fitting engine, replicated to the qpcR package}

\description{
The same function as \code{\link[drc]{drm}}, substitutes the former function \code{multdrc} and has been
 replicated to 'freeze' its working status for qpcR.
}

\usage{
drmfit(formula, curveid, weights, data = NULL, fct, 
       na.action = na.omit, robust = "mean")
}

\arguments{
\item{formula}{a symbolic description of the model to be fit. Either of the form 'response ~ dose' or as a data frame with response values in first column and dose values in second column.}
\item{curveid}{a numeric vector or factor containing the grouping of the data.}
\item{weights}{a numeric vector containing weights.}
\item{data}{an optional data frame containing the variables in the model.}
\item{fct}{one of the available functions of the 'drc' package.}
\item{na.action}{a function which indicates what should happen when the data contain 'NA's.}
\item{robust}{a character string specifying the rho function for robust estimation. See \code{\link[drc]{drm}}.}  
}


\details{
Please use \code{\link{pcrfit}} for easy fitting of single curves
 and \code{drmfit} in combination with \code{\link{repform}} if multiple replicates 
 need to be plotted or analyzed.
}

\value{
An object of class 'drc'.
}

\author{
Christian Ritz and Jens C. Streibig,\cr
modifications by Andrej-Nikolai Spiess
}

\examples{
### on single run data
m <- drmfit(F1.1 ~ Cycles, data = reps, fct = l5())
efficiency(m)
### on replicate data 
repData <- repform(reps[,1:5], c(0,1,1,1,1))
m <- drmfit(values ~ Cycles, curveid = Curve, data = repData, fct = l5())
### plotmean with errorbars
pcrplot(m, type = "errbar")   
}   

\keyword{models}
\keyword{nonlinear}
