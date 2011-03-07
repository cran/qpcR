\name{pcrbatch}
\alias{pcrbatch}

\encoding{latin1}

\title{Batch calculation of qPCR efficiency and other qPCR parameters}

\description{
This function batch calculates the results obtained from \code{\link{efficiency}}, \code{\link{sliwin}} and \code{\link{expfit}} for sigmoidal models or the coefficients (including \code{D0}) from a \code{mak3} model
 on a dataframe containing many qPCR runs.
The input can also be a list obtained from \code{\link{modlist}}, which simplifies things in many cases.
The output is a dataframe with the estimated parameters and model descriptions.
Very easy to use on datasheets containing many qPCR runs, i.e. as can be imported from Excel.  
The result is automatically copied to the clipboard.
}

\usage{
pcrbatch(x, cyc = 1, fluo = NULL, model = l4, remove = FALSE, 
         type = "cpD2", opt = FALSE,  norm = FALSE, backsub = NULL, 
         smooth = c("none", "tukey", "lowess", "supsmu", "spline"), 
         span = 0.1, factor = 1, do.mak = FALSE, opt.method =  "all", 
         nls.method = "all", sig.level = 0.05, crit = "ftest", 
         group = NULL, plot = TRUE, ...) 
}

\arguments{
  \item{x}{a dataframe containing the qPCR raw data from the different runs or a list obtained from \code{\link{modlist}}.}
  \item{cyc}{the column containing the cycle data. Defaults to first column.}
  \item{fluo}{the columns (runs) to be analyzed. If \code{NULL}, all runs will be considered.}
  \item{model}{the model to be used.}
  \item{remove}{logical. If \code{TRUE}, runs that failed to be fitted are not included in the output.}
  \item{type}{the point on the amplification curve from which the efficiency is estimated. See \code{\link{efficiency}}.}
  \item{opt}{logical. Should model optimization take place? If \code{TRUE}, model selection is applied, with criteria as in \code{crit}.}
  \item{norm}{logical. Normalization of the raw data within [0, 1].}
  \item{backsub}{background subtraction. If \code{NULL}, not applied. Otherwise, a numeric sequence such as \code{1:10}. See 'Details'.} 
  \item{smooth}{which curve smoothing method to use. See Details in \code{\link{modlist}}.} 
  \item{span}{the smoothing span for \code{\link{lowess}} or \code{\link{supsmu}}, if selected.} 
  \item{factor}{a constant multiplication factor for the raw qPCR data.}
  \item{do.mak}{logical. If \code{TRUE}, a \code{\link{mak3}} model is fit and the coefficients attached. Time intensive, hence omitted in default.}
  \item{opt.method}{see \code{\link{pcrfit}}.}
  \item{nls.method}{see \code{\link{pcrfit}}.}
  \item{sig.level}{see \code{\link{mselect}}.}     
  \item{crit}{the criterium for model selection. See \code{\link{mselect}}.}
  \item{group}{a vector containing the grouping for possible replicates.} 
  \item{plot}{logical. If \code{TRUE}, the single runs are plotted from the internal 'modlist' for diagnostics.}
  \item{...}{other parameters to be passed to downstream methods.}
}

\details{
The qPCR raw data should be arranged with the cycle numbers in the first column with the name "Cycles".
All subsequent columns must be plain raw data with sensible column descriptions.  
If replicates are defined by \code{group}, the output will contain a numbering of groups (i.e. "group_1" for the first replicate group).
The model selection process is optional, but we advocate using this for obtaining better parameter estimates.
Normalization has been described to improve certain qPCR analyses, but this has still to be independently evaluated. 
Background subtraction is done by averaging the \code{backsub} cycles of the run and subtracting this from all data points.
In case of unsuccessful model fitting, the names are tagged by *NAME* and the column has no values (if \code{remove = FALSE}). However, if \code{remove = TRUE}, the
 failed runs are excluded from the output.
}

\value{
A dataframe with the results in columns containing the calculated values, fit parameters and (tagged) model name together with the different methods used as the name prefix.
 A plot shows a plot matrix of all amplification curves/sigmoidal fits and failed amplifications marked with asterisks.
} 

\seealso{
The function \code{\link{modlist}} for creating a list of models, which is used internally by \code{pcrbatch}.
}

\author{
Andrej-Nikolai Spiess
}

\note{
IMPORTANT: When subsequent use of \code{\link{ratiocalc}} is desired, use \code{pcrbatch} on the single run level with \code{group = NULL} and \code{remove = FALSE}
 (so that \code{\link{ratiocalc}} can automatically delete the failed runs from its \code{group} definition), otherwise error propagation will fail.
}
     
     
\examples{
## complete dataset
\dontrun{
temp <- pcrbatch(reps)
}

## first 4 runs and return parameters of fit
## do background subtraction using the first 8 cycles
res1 <- pcrbatch(reps, fluo = 2:5, backsub = 1:8)

##  first 8 runs, with 4 replicates each, l5 model
res2 <- pcrbatch(reps, fluo = 2:9, model = l5, group = c(1,1,1,1,2,2,2,2))

## using model selection (likelihood ratio) on the first 4 runs, 
## run 1+2 are replicates
res3 <- pcrbatch(reps, fluo = 2:5, group = c(1,1,2,3), opt = TRUE, crit = "ratio")

## converting a 'modlist' to 'pcrbatch'
ml <- modlist(reps, 1, 2:5, b3)
res4 <- pcrbatch(ml)  

\dontrun{
## fitting a 'mak3' mechanistic model
res5 <- pcrbatch(reps, do.mak = TRUE)
View(res5)
}
    
}

\keyword{models}
\keyword{nonlinear}
