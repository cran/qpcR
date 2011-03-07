\name{pcrfit}
\alias{pcrfit}

\encoding{latin1}

\title{Function for qPCR model fitting}

\description{
This is the main workhorse function of the qpcR package that fits one of the available models to qPCR data using nonlinear least-squares fitting from \code{\link{nls}}, 
 with sensible starting parameters obtained from either \code{\link{nls.lm}} (default) or \code{\link{optim}}.
}

\usage{
pcrfit(data, cyc = 1, fluo, model = l4, do.optim = TRUE,  
       opt.method = "all", nls.method = "all", 
       start = NULL, robust = FALSE, weights = NULL, 
       verbose = TRUE, ...)
}

\arguments{
  \item{data}{the name of the dataframe containing the qPCR runs.}
  \item{cyc}{the column containing the cycle data. Defaults to 1.}
  \item{fluo}{the column(s) containing the raw fluorescence data of the run(s). If more than one column is given, the model will be built with the replicates. See 'Details' and 'Examples'.}
  \item{model}{the model to be used for the analysis. Defaults to l4.} 
  \item{do.optim}{if \code{FALSE}, refinement of starting values by any of the optimization methods (see 'Details') will be skipped.}
  \item{opt.method}{all or one of the available optimization methods. See 'Details'.}
  \item{nls.method}{all or one of the available methods in \code{\link{nls}}. See 'Details'.}
  \item{start}{a vector of starting values that can be supplied externally.}
  \item{robust}{logical. If \code{TRUE}, robust nonlinear regression is used. See 'Details'.}   
  \item{weights}{a vector with same length as \code{data} containing possible weights for the nonlinear fit.} 
  \item{verbose}{logical. If \code{TRUE}, fitting and convergence results will be displayed in the console.}
  \item{...}{other parameters to be passed to \code{nls.lm}, \code{optim} or \code{nls}.}     
}

\details{
The fitting procedure works as follows (hopefully ensuring maximum safeness against convergence errors):

1) Approximate starting values are acquired from \code{model$ssfct}.\cr
2) Starting values are refined by any of the methods available in \code{optim} or the Levenberg-Marquardt algorithm (\code{nls.lm}). 
 The default setting is \code{all}, meaning that the optimization is done in the order "LM", "BFGS", "Nelder-Mead", and "SANN".
 If any of the methods successfully converged, the remaining are skipped and the refined starting values transferred to \code{\link{nls}}.  
 The \code{opt.method}s can be combined to tweak the robustness, whereby the starting parameters are passed to
 each succeeding method, i.e. rep("Nelder", 5) will do 5 successive Nelder-Mead optimisations or c("LM", "Nelder") will pass the starting 
 values from "LM" to "Nelder". Levenberg-Marquardt ("LM") is very fast and reliable in many scenarios and is thus the default first method.\cr 
3) The refined starting values from 2) are fitted with \code{\link{nls}} using the order "port", "default" (Gauss-Newton) and "plinear".
 Again, if any of the methods converged, the remaining are skipped.    

This function is to be used at the single run level or on replicates (by giving several columns). The latter will built a model based on the replicate values,
 but only one! If many models should be built on a cohort of replicates, use \code{\link{modlist}} and \code{\link{replist}}.

The output from the optim methods is checked by ensuring all eigenvalues from the hessian are positive, otherwise a notice will occur and the next method in queue is used.  

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
  \item{opt.method}{the successful (convergence) optimization methods.}
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
\emph{SIAM J Appl Math} (1963), \bold{11}: 431-441.\cr
}

\examples{
## simple l4 fit of F1.1 of the 'reps' dataset
m1 <- pcrfit(reps, 1, 2, l4) 
plot(m1)

## using BFGS and Nelder from 'optim'
pcrfit(reps, 1, 2, l5, opt.method = c("BFGS", "Nelder-Mead"))

## skip 'optim' method and supply
## own starting values
pcrfit(reps, 1, 2, l4, do.optim = FALSE, start = c(-5, -0.05, 11, 16)) 

## make a robust model
pcrfit(reps, 1, 2, l4, robust = TRUE)

## make a replicate model
m2 <- pcrfit(reps, 1, 2:5, l5)
plot(m2)

\dontrun{
## fit a mechanistic mak3 model
m3 <- pcrfit(reps, 1, 2, mak3)
plot(m3)
}    
}

\keyword{models}
\keyword{nonlinear}
