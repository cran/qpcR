\name{pcrfit}
\alias{pcrfit}

\encoding{latin1}

\title{Function for qPCR model fitting}

\description{
This is the main workhorse function of the qpcR package that fits one of the available models to qpcR data using non-linear least
 squares fitting from \code{\link{nls}}, with sensible starting parameters obtained from \code{\link{optim}}.
}

\usage{
pcrfit(data, cyc = 1, fluo, model = l4, do.optim = TRUE,  
       opt.method = rep("Nelder", 5), nls.method = "port", 
       start = NULL, ...)
}

\arguments{
  \item{data}{the name of the dataframe containing the qpcR runs.}
  \item{cyc}{the column number with the cycle data. Defaults to 1.}
  \item{fluo}{the column containing the raw fluorescence data of the run.}
  \item{model}{the model to be used for the analysis, defaults to l4.} 
  \item{do.optim}{if \code{FALSE}, refinement of starting values by \code{\link{optim}} will be skipped.}
  \item{opt.method}{one of the available methods in \code{\link{optim}} or "GA". See 'Details'.}
  \item{nls.method}{one of the available methods in \code{\link{nls}}. Default is \code{"port"}, which is quite robust.}
  \item{start}{a vector of starting values that can be supplied externally.}
  \item{...}{other parameters to be passed to \code{optim}, \code{genoud} or \code{nls}.}     
}

\details{
The fitting procedure works as follows (hopefully ensuring maximum robustness against fitting errors):

1) Approximate starting values are acquired from \code{model$ssfct}.\cr
2) Starting values are refined by any of the methods available in \code{optim} or the \code{genoud} method
 from the 'rgenoud' package. The \code{opt.method}s can be combined to tweak the robustness, whereby the starting parameters are passed to
  each succeeding method, i.e. rep("Nelder", 5) will do 5 successive Nelder-Mead optimisations or c("GA", "Nelder") will pass the starting 
   values from "GA" to "Nelder". If problems arise, "GA" has shown to be very robust (but slow!) in the refinement of starting values.
\code{opt.method = c("BFGS", rep("Nelder, 3"))} works also quite well...\cr
3) One of the possible methods from \code{\link{nls}} is then applied with the starting values obtained from 2).    

This function is to be used at the single run level. Otherwise use \code{\link{pcrbatch}} or \code{\link{modlist}}.

The output from the optim methods is checked by ensuring all eigenvalues from the hessian are positive, otherwise a notice will occur.  
}

\value{
A model of class 'nls' and 'pcrfit' with the following items attached:\cr
  \item{DATA}{the initial data used for fitting.}
  \item{MODEL}{the model used for fitting.}  
  \item{parMat}{the trace of the starting values for each applied method. Can be used to track problems.}
  \item{opt.method}{the parameter \code{opt.method}.}
} 

\author{
Andrej-Nikolai Spiess
}


\references{
Ritz C & Streibig JC. Bioassay analysis using R. \emph{J Stat Soft} (2005), \bold{12}: 1-22.
}

\examples{
## simple l4 fit of F1.1 of the 'reps' dataset
pcrfit(reps, 1, 2, l4) 

## same with five-parameter model
## use "GA" method for optim
pcrfit(reps, 1, 2, l5, opt.method = "GA")

## skip 'optim' method and supply
## own starting values
pcrfit(reps, 1, 2, l4, do.optim = FALSE, start = c(-5, -0.05, 11, 16)) 
}

\keyword{models}
\keyword{nonlinear}
