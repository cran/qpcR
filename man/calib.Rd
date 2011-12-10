\name{calib}
\alias{calib}

\title{Calculation of qPCR efficiency by dilution curve analysis}

\description{
This function calculates the PCR efficiency from a classical qPCR dilution experiment. The threshold cycles are plotted against the logarithmized concentration (or dilution) values, a linear regression line is fit and the efficiency calculated by \eqn{E = 10^{\frac{-1}{slope}}}. A graph is displayed with the raw values plotted with the threshold cycle and the linear regression curve.  
The threshold cycles are calculated either by some arbitrary fluorescence value (i.e. as given by the qPCR software) or calculated from the second derivative maximum of the first (highest concentration) curve. If values to be predicted are given, they are calculated from the curve and also displayed within. An iterative search of the optimal threshold border can be conducted, thereby going through all slopes and intercepts and selecting the combination that minimizes the AIC of the acquired linear regression curves. See 'Details' for more information on this (maybe controversial) procedure.   
} 

\usage{
calib(refcurve, predcurve = NULL, thresh = "cpD2", term = NULL, 
      dil = NULL, plot = TRUE, plot.map = TRUE,
      conf = 0.95, opt = c("none", "inter", "slope"),       
      opt.step = c(50, 50), stop.crit = c("outlier", "midpoint"), 
      quan = 0.5, slope = NULL, count = 1)
}

\arguments{
 \item{refcurve}{a 'modlist' containing the curves for calibration.}
 \item{predcurve}{an (optional) 'modlist' containing the curves for prediction.}
 \item{thresh}{the fluorescence value from which the threshold cycles are defined. Either "cpD2" or a numeric value.} 
 \item{term}{an (optional) numeric value for the terminating intercept. See 'Details'.}
 \item{dil}{a vector with the concentration (or dilution) values corresponding to the calibration curves.}
 \item{plot}{logical. Should the optimization process be displayed? If \code{FALSE}, only values are returned.}
 \item{plot.map}{logical. Should a final heatmap display from the goodness-of-fit of all iterations be displayed?}
 \item{conf}{the confidence interval. Defaults to 95\%, can be omitted with \code{NULL}.}
 \item{opt}{type of optimization. See 'Details'.}
 \item{opt.step}{a two-element vector. Number of iterations for the intercept as first item, number of slope iterations as second.}
 \item{stop.crit}{the stopping criterium for the iterations. Default is the first outlier cycle of the highest dilution curve.}
 \item{quan}{the top quantile of iterations to be shown in the heatmap.}
 \item{slope}{a slope to be defined for the threshold line. Mostly used internally by the iteration process.}
 \item{count}{internal counter for recursive purposes. Not to be altered.}  
}

\details{                    
The iterative function conducts a search through all combinations of slope and intercept. For each iteration, either \eqn{R^2} or AIC of the resulting calibration curves are collected, and finally the combination is selected that minimized the AIC. The function goes through all combinations as to avoid local maxima that are likely to happen in this approach. The different settings for \code{opt} are:

"none"  only second derivative maximum or single threshold value.\cr
"inter" iterate along the y-axis intercept.\cr
"slope" iterate along y-axis intercept and slope values.

The paradigm is such that the iterations will start at the second derivative of the first (lowest dilution; highest copy number) curve and
 terminate at the outlier cycle of the last (highest dilution; lowest copy number) curve. Alternatively a y-value can be
 given with \code{term} for the termination threshold. The number of iterations can be defined by \code{opt.step} but the default values
 usually give reasonable results. Not to forget, an iterative search only throughout all intercepts can be chosen, as well as a classical
 approach using the second derivative maximum of the first curve or a defined threshold value from the qPCR software, to be defined in 
 \code{thresh}. See 'Examples'.
}

\value{
  A list with the following components:
  \item{refcyc}{the calculated threshold cycles for the calibration curves.}
  \item{refcyc.conf}{the confidence intervals for \code{refcyc}.}
  \item{predcyc}{the calculated threshold cycles for the prediction curves.}  
  \item{predconc}{the predicted concentrations.}
  \item{predconc.conf}{the confidence intervals for \code{predconc}.}    
  \item{eff}{the efficiency as calculated from the calibration curve.}
  \item{aic}{the AIC value of the linear fit.}
  \item{aicc}{the corrected AIC value of the linear fit.}   
  \item{rsq}{The \eqn{R^2} of the linear fit.} 
  \item{rsq.ad}{The adjusted \eqn{R^2_{adj}} of the linear fit.}   
  \item{aicMat}{a matrix with the calibration curve AIC of each iteration.}  
  
ATTENTION: If iterations were used, the values reflect the analysis of the best fit! 
}

\author{
Andrej-Nikolai Spiess
}       

\references{
A new mathematical model for relative quantification in real-time RT-PCR.\cr
Pfaffl MW.\cr
\emph{Nucleic Acids Res} (2001), \bold{29}:e45.
}

\examples{
## Define calibration curves,
## dilutions (or copy numbers) 
## and curves to be predicted.
## Do background subtraction using
## average of first 8 cycles
CAL <- modlist(reps, fluo = c(2, 6, 10, 14, 18, 22), baseline = 1:8)
COPIES <- c(100000, 10000, 1000, 100, 10, 1)
PRED <- modlist(reps, fluo = c(3, 7, 11), baseline = 1:8)

## conduct normal quantification using
## the second derivative maximum of 
## first curve
res1 <- calib(refcurve = CAL, predcurve = PRED, thresh = "cpD2", dil = COPIES) 

## using a defined treshold value
res2 <- calib(refcurve = CAL, predcurve = PRED, thresh = 0.5, dil = COPIES) 

\dontrun{
## using replicates for reference curve
ml1 <- modlist(reps, model = l4)
DIL <- rep(10^(6:0), each = 4)
res1 <- calib(refcurve = ml1, dil = DIL)

## iterating the intercept with 50 steps
par(ask = FALSE)
res2 <- calib(refcurve = CAL, predcurve = PRED, dil = COPIES, opt = "inter",
              opt.step = c(50, 0))

## iterating the intercept/slope with 20 steps
par(ask = FALSE)
res3 <- calib(refcurve = CAL, predcurve = PRED, dil = COPIES, opt = "slope", 
              opt.step = c(20, 20)) 
}
}

\keyword{models}
\keyword{nonlinear}
