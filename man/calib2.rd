\name{calib2}
\alias{calib2}

\title{Calculation of qPCR efficiency by dilution curve analysis and an iterative search of the optimal threshold border}

\description{
\code{calib2} is an advanced version of \code{\link{calib}}! An iterative search of the optimal threshold border is conducted,
 thereby going through all slopes and intercepts and selecting the combination that minmimizes the AIC of the acquired
 linear regression curves. See 'Details' for more information on this (maybe controversial) procedure. 
 Also different to \code{\link{calib}}, this function uses two 'modlist's, one for calibration and one for prediction.  
} 

\usage{
  calib2(refcurve, predcurve = NULL, thresh = "cpD2", term = NULL, 
        dil = NULL, fct = l5(), plot = TRUE, plot.map = TRUE, 
        conf = 0.95, opt = c("none", "inter", "slope"),
        opt.step = c(50, 50), quan = 0.5, slope = NULL, count = 1)
}

\arguments{
 \item{refcurve}{a 'modlist' containing the curves for calibration.}
 \item{predcurve}{an (optional) 'modlist' containing the curves for prediction.}
 \item{thresh}{the fluorescence value from which the threshold cycles are defined. Either "cpD2" or a numeric value.} 
 \item{term}{an (optional) numeric value for the terminating intercept. See 'Details'.}
 \item{dil}{a vector with the concentration (or dilution) values corresponding to the calibration curves.}
 \item{fct}{the model used for the threshold cycle estimation. Any sigmoidal model in the 'drc' package. Defaults to 'l5()'.} 
 \item{plot}{logical. Should the optimization process be displayed? If \code{FALSE}, only values are returned.}
 \item{plot.map}{logical. Should a final heatmap display from the goodness-of-fit of all iterations be displayed?}
 \item{conf}{the p-value for the confidence interval. Defaults to 0.05, can be omitted with \code{NULL}.}
 \item{opt}{type of optimization. See 'Details'.}
 \item{opt.step}{a two-element vector. Number of iterations for the intercept as first item, number for slope iterations as second.}
 \item{quan}{the top quantile of iterations to be shown in the heatmap.}
 \item{slope}{a slope to be defined for the threshold line. Mostly used internally by the iteration process.}
 \item{count}{internal counter for recursive purposes. Not to be altered.}  
}

\details{
This function conducts an iterative search through all combinations of slope and intercept. For each iteration, either R-square or AIC
 of the resulting calibration curves are collected, and finally the combination is selected that minimized the AIC.
The function goes through all combinations as to avoid local maxima that are likely to happen in this approach.
The different settings for \code{opt} are:

"none"  only second derivative maximum or single threshold value.\cr
"inter" iterate along the y-axis intercept.\cr
"slope" iterate along y-axis intercept and slope values.

The paradigm is such that the iterations will start at the second derivative of the first (lowest dilution; highest copy number) curve and
 terminate at the outlier cycle of the last (highest dilution; lowest copy number) curve. Alternatively a defined y-value can be
 defined with \code{term} for the termination threshold. The number of iterations can be defined by \code{opt.step} but the default values
  usually give reasonable results. Not to forget, an iterative search ONLY throughout all intercepts can be chosen, as well as a classical
  approach using only the second derivative maximum of the first curve or a defined threshold value from the qPCR software, to be defined in 
   \code{thresh}. See 'Examples'.
}

\value{
  A list with the following components:
  \item{ref.cyc}{the calculated threshold cycles for the calibration curves.}
  \item{pred.cyc}{the calculated threshold cycles for the curves to be predicted.}
  \item{pred.logconc}{the predicted concentrations on a log scale.}
  \item{pred.conc}{the predicted concentrations on a linear scale.}
  \item{ref.conf}{the confidence values for the calibration curve points. In the same order as \code{ref.cyc}.}
  \item{pred.logconf}{the confidence values for the prediction curve points, in the same order as \code{pred.cyc} and on a log scale.}
  \item{pred.conf}{the confidence values for the prediction curve points, in the same order as \code{pred.cyc} and on a linear scale.}
  \item{eff}{the efficiency as calculated from the calibration curve.}
  \item{aic}{the AIC value of the linear fit.}
  \item{rsq}{The r-square of the linear fit.} 
  \item{aicMat}{a matrix with the calibration AIC of each iteration.}  
  
  ATTENTION: If iterations were used, the values reflect the analysis of the best fit! 
}

\author{
  Andrej-Nikolai Spiess
}

\references{
A paper describing the benefits of this approach will follow...
}

\examples{
## Define calibration curves,
## dilutions (or copy numbers) 
## and curves to be predicted.
## Do background subtraction using
## average of first 8 cycles
CAL <- modlist(reps, c(2, 6, 10, 14, 18, 22), backsub = 1:8)
COPIES <- c(100000, 10000, 1000, 100, 10, 1)
PRED <- modlist(reps, c(3, 7, 11), backsub = 1:8)
## conduct normal quantification using
## the second derivative maximum of 
## first curve
res <- calib2(refcurve = CAL, predcurve = PRED, thresh = "cpD2", dil = COPIES) 
## using a defined treshold value
res <- calib2(refcurve = CAL, predcurve = PRED, thresh = 0.5, dil = COPIES) 
## iterating only the intercept with 20 steps
res <- calib2(refcurve = CAL, predcurve = PRED, dil = COPIES, opt = "inter", 
              opt.step = c(20, 0)) 
## iterating only the intercept/slope with 10 steps
res <- calib2(refcurve = CAL, predcurve = PRED, dil = COPIES, opt = "slope", 
              opt.step = c(5, 10)) 
}

\keyword{models}
\keyword{nonlinear}
