\name{pcrbatch}
\alias{pcrbatch}

\encoding{latin1}

\title{Batch calculation of qPCR efficiency and several other important qPCR parameters with different methods}

\description{
This function batch calculates the results obtained from \code{\link{efficiency}}, \code{\link{sliwin}}, and \code{\link{expfit}} on a dataframe containing many qPCR runs.
The output is a dataframe with the estimated parameters and model description.
Very easy to use on datasheets containing many qPCR runs, i.e. as can be imported from Excel.  
The data is automatically copied to the clipboard.
}

\usage{
pcrbatch(x, cols = 2:ncol(x), group = NULL, model = l4(), type = "cpD2", 
	 opt = FALSE, smooth = c("none", "tukey", "lowess"), norm = FALSE, fact = 1, 
	 ave = c("mean", "median"), plot = FALSE, retPar = FALSE, ...) 
}

\arguments{
  \item{x}{a dataframe containing the qPCR raw data from the different runs.}
  \item{cols}{the columns (runs) to be analyzed. Defaults to every run in the dataframe.}
  \item{group}{a vector containing the grouping for possible replicates.} 
  \item{model}{the model to be used for the analysis, with the usual 'drc' package nomenclature.}
  \item{type}{the point on the amplification curve from which the efficiency is estimated. See \link{efficiency}.}
  \item{opt}{logical. Should model optimization take place? If \code{TRUE}, model selection is applied.}
  \item{smooth}{the smoothing algorithm for the data. Defaults to Tukey´s running median. Non-parametric lowess smoothing with "lowess" is another option.}
  \item{norm}{logical. Normalization of the raw data such that the highest value (plateau phase) is 1. See references.}
  \item{fact}{a constant multiplication factor for the raw qPCR data.}
  \item{ave}{averaging method for replicates. Defaults to "mean", another option is "median".}
  \item{plot}{logical. If \code{TRUE}, a graphical analysis of the qPCR data is displayed, otherwise (default) a simple console output is displayed.}
  \item{retPar}{logical. Should the parameters from the fit be included in the output?}
  \item{...}{other parameters to be passed to \code{\link{mchoice}}.}
}

\details{
The qPCR raw data should be arranged with the cycle numbers in the first column with the name "Cycles".
All subsequent columns must be plain raw data with sensible column descriptions.
If replicates are defined, the output will contain a numbering of groups (i.e. "group1" for the first replicate group).
The model selection process is optional, but we advocate using this for obtaining better parameter estimates.
Normalization has been described to improve certain qPCR analyses, but this has still to be independently evaluated. 
}

\value{
  A dataframe with the results in columns containing the calculated values with descriptions and the method used as the name prefix.
} 

\author{
  Andrej-Nikolai Spiess
}

\note{
When subsequent use of \code{\link{ratiocalc}} is desired, use \code{pcrbatch} on the single run level with \code{group = NULL}.
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
temp <- pcrbatch(reps, 2:4, retPar = TRUE)

##  first 8 runs, with 4 replicates each
temp <- pcrbatch(reps, 2:9, c(1,1,1,1,2,2,2,2))

## using model selection (likelihood ratio) on the first 4 runs, run 1+2 are replicates
temp <- pcrbatch(reps, 2:5, c(1,1,2,3), opt = TRUE, crit = "ratio")
}

\keyword{models}
\keyword{nonlinear}
