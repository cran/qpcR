\name{pcrsim}
\alias{pcrsim}

\title{Simulation of sigmoidal qPCR data with goodness-of-fit analysis for different models}

\description{Simulated sigmoidal qPCR curves are generated from an initial model to which some
 user-defined noise is added. One or more models can then be fit to this random data and goodness-of-fit (GOF)
 measures are calculated for each of the models. This is essentially a Monte-Carlo approach testing for
 the best model in dependence to some noise structure in sigmodal models.
}

\usage{
pcrsim(cyc = 1:30, model = l4, par = NULL, nsim = 10,        
       error = 0.02, errfun = function(y) 1, plot = TRUE,
       fitmodel = NULL, select = FALSE, 
       statfun = function(y) mean(y, na.rm = TRUE), ...)
}

\arguments{
  \item{cyc}{the number of cycles to be simulated.}
  \item{model}{the initial model used for generating the simulated curves.} 
  \item{par}{a numeric vector (required) with the parameter values for \code{model}.} 
  \item{nsim}{the number of simulated curves.}
  \item{error}{the gaussian error used for the simulation. See 'Details'.}
  \item{errfun}{an optional function for the error distribution. See 'Details'.}
  \item{plot}{should the simulated and fitted curves be displayed?}
  \item{fitmodel}{a model or model list to test against the initial model.} 
  \item{select}{if \code{TRUE}, a matrix is returned with the best model in respect to each of the GOF measures.}  
  \item{statfun}{a function to be finally applied to all collected GOF measures, default is the average.} 
  \item{...}{other parameters to be passed on to \code{\link{plot}} or \code{\link{pcrfit}}.}
}

\details{
The value defined under \code{error} is just the standard deviation added plainly to each y_i value from the initial model, thus generating
 a dataset with random homoscedastic noise. With aid of \code{errfun}, the distribution of the error along the y_i values can
 be altered and therefore be used to generate heteroscedastic variance along the curve, so that the standard deviation is a function
 of the magnitude.
 
Example:\cr
\code{errfun = function(y) 1}\cr
same variance for all y_i, as is.\cr

\code{errfun = function(y) y}\cr
variance as a function of the y-magnitude.\cr

\code{errfun = function(y) 1/y}\cr
variance as an inverse function of the y-magnitude.

For the effect, see 'Examples'.
}

\value{
A list containing the following items:
  \item{cyc}{same as in 'arguments'.}
  \item{fluoMat}{a matrix with the simulated qPCR data in columns.}
  \item{coefList}{a list with the coefficients from the fits for each model, as subitems.}
  \item{gofList}{a list with the GOF measures for each model, as subitems.} 
  \item{statList}{a list with the GOF measures summarized by \code{statfun} for each model, as subitems.} 
  \item{modelMat}{if \code{select = TRUE}, a matrix with the best model for each GOF measure and each simulation.}   
}

\author{
Andrej-Nikolai Spiess
}  

\examples{
## generate initial model
m <- pcrfit(reps, 1, 2, l4)

## simulate homoscedastic error
## and test initial model, w3 and l5 
## model on data
res <- pcrsim(cyc = 1:30, model = l4, par = coef(m),
              error = 0.2, nsim = 20, fitmodel = list(l4, w3, l5))

## use heteroscedastic noise typical for 
## qPCR: more noise at lower fluorescence
\dontrun{
res2 <- pcrsim(cyc = 1:30, model = l4, par = coef(m),
              error = 0.01, errfun = function(y) 1/y,
              nsim = 20, fitmodel = list(l4, w3, l5))
}

## get 95\% confidence interval for 
## the models GOF in question (l4, w3, l5) 
\dontrun{
res <- pcrsim(cyc = 1:30, model = l4, par = coef(m),
              error = 0.2, nsim = 20, fitmodel = list(l4, w3, l5),
              statfun = function(y) quantile(y, c(0.025, 0.975)))
res$statList  
}  

## show best model for each simulation
## based on different GOF measures
\dontrun{
m2 <- pcrfit(reps, 1, 2, l3)
res <- pcrsim(cyc = 1:30, model = l3, par = coef(m2),
              error = 0.2, nsim = 200, fitmodel = list(l3, l4, l5),
              select = TRUE)
res$modelMat
} 
}  





\keyword{models}
\keyword{nonlinear}
