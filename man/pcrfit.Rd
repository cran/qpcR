\name{pcrfit}
\alias{pcrfit}

\encoding{latin1}

\title{Function for qPCR model fitting}

\description{
This is the main workhorse function of the qpcR package that fits one of the available models to qPCR data using nonlinear least-squares fitting from \code{\link{nls}}, 
 with sensible starting parameters obtained from either \code{nls.lm}, \code{\link{optim}} or \code{genoud}.
}

\usage{
pcrfit(data, cyc = 1, fluo, model = l4, do.optim = TRUE,  
       opt.method = "LM", nls.method = "port", 
       start = NULL, robust = FALSE, control = nls.control(), 
       weights, ...)
}

\arguments{
  \item{data}{the name of the dataframe containing the qPCR runs.}
  \item{cyc}{the column containing the cycle data. Defaults to 1.}
  \item{fluo}{the column containing the raw fluorescence data of the run.}
  \item{model}{the model to be used for the analysis. Defaults to l4.} 
  \item{do.optim}{if \code{FALSE}, refinement of starting values by \code{nls.lm}, \code{optim} or \code{nls.lm} will be skipped.}
  \item{opt.method}{one of the available refinement methods. See 'Details'.}
  \item{nls.method}{one of the available methods in \code{\link{nls}}. Default is \code{"port"}, which works quite well.}
  \item{start}{a vector of starting values that can be supplied externally.}
  \item{robust}{logical. If \code{TRUE}, robust nonlinear regression is used. See 'Details'.}
  \item{control}{an optional list of control settings for \code{\link{nls}} or \code{qpcR:::rnls}.}
  \item{weights}{a vector of weights for nonlinear fitting. Must be same lengths as the data.}
  \item{...}{other parameters to be passed to \code{optim}, \code{genoud} or \code{nls}.}     
}

\details{
The fitting procedure works as follows (hopefully ensuring maximum safeness against convergence errors):

1) Approximate starting values are acquired from \code{model$ssfct}.\cr
2) Starting values are refined by any of the methods available in \code{optim}, the \code{genoud} method
 from the 'rgenoud' package or the Levenberg-Marquardt algorithm (\code{nls.lm}). The \code{opt.method}s can be combined to tweak the robustness, whereby the starting parameters are passed to
  each succeeding method, i.e. rep("Nelder", 5) will do 5 successive Nelder-Mead optimisations or c("GA", "Nelder") will pass the starting 
   values from "GA" to "Nelder". If problems arise, "GA" has shown to be very robust (but slow!) in the refinement of starting values.
 Levenberg-Marquardt ("LM") is very fast and relatively reliable in many scenarios and is thus the default.\cr 
3) One of the possible methods from \code{\link{nls}} is then applied with the starting values obtained from 2).    

This function is to be used at the single run level. Otherwise use \code{\link{pcrbatch}} or \code{\link{modlist}}.

The output from the optim methods is checked by ensuring all eigenvalues from the hessian are positive, otherwise a notice will occur.  

If \code{robust = TRUE}, robust nonlinear fitting will be used. To do this, the internal function \code{qpcR:::rnls} is called
 which is a modification of the \code{nlrob} function of the 'robustbase' package. Modifications were done such that all available generic functions
  for objects of class 'nls' can be used on the output of \code{qpcR:::rnls}, such as \code{predict}, \code{confint} etc.
}

\value{
A model of class 'nls' and 'pcrfit' with the following items attached:\cr
  \item{DATA}{the initial data used for fitting.}
  \item{MODEL}{the model used for fitting.} 
  \item{call2}{the call to \code{pcrfit}.} 
  \item{parMat}{the trace of the starting values for each applied method. Can be used to track problems.}
  \item{opt.method}{the parameter \code{opt.method}.}
} 

\author{
Andrej-Nikolai Spiess
}


\references{
Bioassay analysis using R.\cr
Ritz C & Streibig JC.\cr
\emph{J Stat Soft} (2005), \bold{12}: 1-22.\cr

A Method for the Solution of Certain Problems in Least Squares.\cr
K. Levenberg.\cr
\emph{Quart Appl Math} (1944), \bold{2}: 164-168.\cr

An Algorithm for Least-Squares Estimation of Nonlinear Parameters.\cr
D. Marquardt.\cr
\emph{SIAM J Appl Math} (19639), \bold{11}: 431-441.
}

\examples{
## simple l4 fit of F1.1 of the 'reps' dataset
pcrfit(reps, 1, 2, l4) 

## same with five-parameter model
## use "GA" method for optim
pcrfit(reps, 1, 2, l5, opt.method = "GA")

## using BFGS and Nelder from 'optim'
pcrfit(reps, 1, 2, l5, opt.method = c("BFGS", "Nelder"))

## skip 'optim' method and supply
## own starting values
pcrfit(reps, 1, 2, l4, do.optim = FALSE, start = c(-5, -0.05, 11, 16)) 

## make a robust model
pcrfit(reps, 1, 2, l4, robust = TRUE) 

## use weights
pcrfit(reps, 1, 2, l4, weights = 49:1)  
}

\keyword{models}
\keyword{nonlinear}
