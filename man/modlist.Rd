\name{modlist}
\alias{modlist}

\title{Create nonlinear models from a dataframe and coerce them into a list}

\description{
Essential function to create a list of nonlinear models from the columns (runs) of a qPCR dataframe. This function houses different methods for curve transformation prior to fitting, such as normalization in [0, 1], smoothing, baseline subtraction etc. Runs that failed to fit or that have been identified as kinetic outliers (by default: lack of sigmoidal structure) can be removed automatically as well as their entries in an optionally supplied label vector. 
}

\usage{
modlist(x, cyc = 1, fluo = NULL, model = l4, check = "uni2", 
        checkPAR = parKOD(), remove = c("none", "fit", "KOD"), 
        exclude = NULL, labels = NULL, norm = FALSE, baseline = FALSE,
        basefac = 1, smooth = c("none", "smooth", "lowess", "supsmu", "spline"), 
        smoothPAR = list(span = 0.1), factor = 1, opt = FALSE,
        optPAR = list(sig.level = 0.05, crit = "ftest"), verbose = TRUE, ...)
}

\arguments{
  \item{x}{a dataframe containing the qPCR data or a single qPCR run of class 'pcrfit'.}
  \item{cyc}{the column containing the cycle data. Defaults to first column.}
  \item{fluo}{the column(s) (runs) to be analyzed. If \code{NULL}, all runs will be considered.}
  \item{model}{the model to be used for all runs.} 
  \item{check}{the method for kinetic outlier detection. Default is check for sigmoidal structure, see \code{\link{KOD}}. To turn off, use \code{NULL}.}
  \item{checkPAR}{parameters to be supplied to the \code{check} method, see \code{\link{KOD}}.}
  \item{remove}{which runs to remove. Either \code{"none"}, those which failed to \code{"fit"} or from the \code{"KOD"} outlier method.}
  \item{exclude}{either \code{""} for samples with missing column names or a regular expression defining columns (samples) to be excluded from \code{modlist}. See 'Details'.} 
  \item{labels}{a vector containing labels, i.e. for defining replicate groups prior to \code{\link{ratiobatch}}.}
  \item{norm}{logical. Should the raw data be normalized within [0, 1] before model fitting?}
  \item{baseline}{type of baseline subtraction. Same as in \code{\link{efficiency}}.}
  \item{basefac}{a factor when using averaged baseline cycles, such as \code{0.95}.}
  \item{smooth}{which curve smoothing method to use. See 'Details'.}
  \item{smoothPAR}{parameters to be supplied to the smoothing functions. See 'Details'.}
  \item{factor}{a multiplication factor for the fluorescence response values (barely useful, but who knows...).}
  \item{opt}{logical. Should model selection be applied to each model?}
  \item{optPAR}{parameters to be supplied to \code{\link{mselect}}.}
  \item{verbose}{logical. If \code{TRUE}, fitting and tagging results will be displayed in the console.}
  \item{...}{other parameters to be passed to \code{\link{pcrfit}}.}
}

\value{
A list with each item containing the model from each column. A \code{names} item (which is tagged by *NAME*, if fitting failed) containing the column name is attached to each model as well as an item \code{isFitted} with either \code{TRUE} (fitting converged) or \code{FALSE} (a fitting error occured). This information is useful when \code{\link{ratiocalc}} is to be applied and unsuccessful fits should automatically removed from the given \code{group} definition. If kinetic outlier detection is selected, an item \code{isOutlier} is attached, defining the run as an outlier (\code{TRUE}) or not (\code{FALSE}). 
}

\details{
The following smoothing methods are available for the fluorescence values:\cr
\code{\link{smooth}}, \code{\link{lowess}}, \code{\link{supsmu}} and \code{\link{spline}}. See documentation there. The author of this package favorizes \code{"supsmu"} with \code{span = 0.1}. In case of unsuccessful model fitting and if \code{remove = "none"} (default), the original data is included in the output, albeit with no fitting information. This is useful since using \code{plot.pcrfit} on the 'modlist' shows the non-fitted runs. If \code{remove = "fit"}, the non-fitted runs are automatically removed and will thus not be displayed. If \code{remove = "KOD"}, by default all runs without sigmoidal structure are removed likewise. If a \code{labels} vector \code{lab} is supplied, the labels from the failed fits are removed and a new label vector \code{lab_mod} is written to the global environment. This way, an initial labeling vector for all samples can be supplied, bad runs and their labels automatically removed and these transfered to downstream analysis (i.e. to \code{\link{ratiobatch}}) without giving errors. \code{exclude} offers an option to exclude samples from the modlist by some regular expression or by using \code{""} for samples with empty column names. See 'Examples'. 
}

\seealso{
\code{\link{pcrbatch}} for batch analysis using different methods.
}

\author{
Andrej-Nikolai Spiess
}

\examples{
## calculate efficiencies and ct values 
## for each run in the 'reps' data, 
## subtract baseline using the 
## first 8 cycles.
ml1 <- modlist(reps, model = l5, baseline = 1:8)
getPar(ml1, type = "curve")

## 'crossing points' for the first 3 runs (normalized)
##  and using best model from Akaike weights
ml2 <- modlist(reps, 1, 2:5, model = l5, norm = TRUE, 
               opt = TRUE, optPAR = list(crit = "weights"))
sapply(ml2, function(x) efficiency(x, plot = FALSE)$cpD2)

## convert a single run to a 'modlist'
m <- pcrfit(reps, 1, 2, l5)
ml3 <- modlist(m)

## using the 'testdat' set
## include failed fits
ml4 <- modlist(testdat, 1, 2:9,  model = l5)
plot(ml4, which = "single")

## remove failed fits and update a label vector
GROUP <- c("g1s1", "g1s2", "g1s3", "g1s4", "g1c1", "g1c2", "g1c3", "g1c4") 
ml5 <- modlist(testdat, 1, 2:9,  model = l5, labels = GROUP, remove = "KOD")
plot(ml5, which = "single")
GROUP_mod

\dontrun{
## use one of the mechanistic models
## get D0 values
ml6 <- modlist(reps, model = mak3)
sapply(ml6, function(x) coef(x)[1])

## exclude first sample in each 
## replicate group of dataset 'reps'
ml7 <- modlist(reps, exclude = ".1")
plot(ml7, which = "single")

## using weighted fitting:
## weighted by inverse residuals
ml8 <- modlist(reps, weights = "1/abs(resid)")
plot(ml8, which = "single")
}
}

\keyword{models}
\keyword{nonlinear}
