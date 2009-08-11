\name{pcrbatch}
\alias{pcrbatch}

\encoding{latin1}

\title{Batch calculation of qPCR efficiency and several other important qPCR parameters with different methods}

\description{
This function batch calculates the results obtained from \code{\link{efficiency}}, \code{\link{sliwin}}, and \code{\link{expfit}} on a dataframe containing many qPCR runs.
The input can also be a list obtained from \code{\link{modlist}}, which simplifies things in many cases.
The output is a dataframe with the estimated parameters and model description.
Very easy to use on datasheets containing many qPCR runs, i.e. as can be imported from Excel.  
The data is automatically copied to the clipboard.
}

\usage{
pcrbatch(x, cols = NULL, model = l4, group = NULL, type = "cpD2",
         opt = FALSE, smooth = c("none", "tukey", "lowess"), norm = FALSE, 
         fact = 1, ave = c("mean", "median"), backsub = NULL, 
         retPar = FALSE, crit, ...) 
}

\arguments{
  \item{x}{a dataframe containing the qPCR raw data from the different runs or a list obtained from \code{\link{modlist}}.}
  \item{cols}{the columns (runs) to be analyzed. If \code{NULL}, all runs will be considered.}
  \item{group}{a vector containing the grouping for possible replicates.} 
  \item{model}{the model to be used.}
  \item{type}{the point on the amplification curve from which the efficiency is estimated. See \link{efficiency}.}
  \item{opt}{logical. Should model optimization take place? If \code{TRUE}, model selection is applied.}
  \item{smooth}{the smoothing algorithm for the data. Either Tukey´s running median or non-parametric lowess smoothing.}
  \item{norm}{logical. Normalization of the raw data within [0, 1]. See references.}
  \item{fact}{a constant multiplication factor for the raw qPCR data.}
  \item{ave}{averaging method for replicates. Defaults to "mean", another option is "median".}
  \item{backsub}{background subtraction. If \code{NULL}, not applied. Otherwise, a numeric sequence such as \code{1:10}. See 'Details'.} 
  \item{retPar}{logical. Should the parameters from the fit be included in the output?}
  \item{crit}{the criterium for model selection. See \code{\link{mselect}}.}
  \item{...}{other parameters to be passed to downstream methods.}
}

\details{
The qPCR raw data should be arranged with the cycle numbers in the first column with the name "Cycles".
All subsequent columns must be plain raw data with sensible column descriptions.
If data of class \code{modlist} is given, all previous data manipulations (i.e. smoothing or background subtraction) are eliminated as the raw data is taken.
Thus, if \code{opt} was applied with \code{modlist}, this has to be repeated by setting \code{opt = TRUE}.
If replicates are defined, the output will contain a numbering of groups (i.e. "group1" for the first replicate group).
The model selection process is optional, but we advocate using this for obtaining better parameter estimates.
Normalization has been described to improve certain qPCR analyses, but this has still to be independently evaluated. 
Background subtraction is done by averaging the \code{backsub} cycles of the run and subtracting this from all data points.
}

\value{
A dataframe with the results in columns containing the calculated values with descriptions and the method used as the name prefix.
} 

\author{
Andrej-Nikolai Spiess
}

\note{
When subsequent use of \code{\link{ratiocalc}} is desired, use \code{pcrbatch} on the single run level with \code{group = NULL},
 otherwise error propagation will fail.
}

\references{
A standard curve based method for relative real time PCR data processing.
Larionov et al., \emph{BMC Bioinformatics}, \bold{6}: 62.
}


\examples{
## complete dataset
\dontrun{
temp <- pcrbatch(reps)
}

## first 4 runs and return parameters of fit
## do background subtraction using the first 8 cycles
temp <- pcrbatch(reps, 2:4, retPar = TRUE, backsub = 1:8)

##  first 8 runs, with 4 replicates each, l5 model
temp <- pcrbatch(reps, 2:9, model = l5, c(1,1,1,1,2,2,2,2))

## using model selection (likelihood ratio) on the first 4 runs, 
## run 1+2 are replicates
temp <- pcrbatch(reps, 2:5, group = c(1,1,2,3), opt = TRUE, crit = "ratio")

## converting a 'modlist' to 'pcrbatch'
ml <- modlist(reps, 1, 2:5, l4)
temp <- pcrbatch(ml)
}

\keyword{models}
\keyword{nonlinear}
