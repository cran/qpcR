\name{calib}
\alias{calib}

\title{Calculation of qPCR efficiency by dilution curve analysis}

\description{
  This function calculates the PCR efficiency from a classical qPCR dilution experiment. The threshold cycles are plotted against the 
  logarithmized concentration (or dilution) values, a linear regression line is fit and the efficiency calculated by \eqn{E = 10^{\frac{-1}{slope}}}.
  A graph is displayed with the raw values plotted with the threshold cycle and the linear regression curve.  
  The threshold cycles are calculated either by some arbitrary fluorescence value (i.e. as given by the qPCR software)
  or calculated from the second derivative maximum of the first (highest concentration) curve. 
  If values to be predicted are given, they are calculated from the curve and also displayed within.
  The confidence interval (p = 0.05) for the regression line is displayed at default but can be omitted. 
  An optimization procedure is implemented which starts from the first derivative maximum and stops at a user defined
  threshold value. For 1000 threshold cycles within this region the Akaike Information Criterion (AIC) is calculated for 
  the regression line of the dilution steps. This method aims to find the threshold cycle with the optimal fit and
  must be tweaked by the user (otherwise fitting errors will occur!).
} 

\usage{
  calib(data, cyc.col = 1, f.col = 2:ncol(data), groups = NULL, thresh = 1, 
        method = "var", dil = NULL, pred = NULL, fct = l5(), plot = TRUE, 
        conf.lev = 0.95, opt.min = NULL, opt.plot = FALSE)
}

\arguments{
 \item{data}{a dataframe containing the qPCR data.}
 \item{cyc.col}{the column containing the cycle numbers. Defaults to the first column.}
 \item{f.col}{the columns containing the raw fluorescence qPCR data.} 
 \item{groups}{a factor defining replicates of the qPCR runs.}
 \item{thresh}{the fluorescence value from which the threshold cycles are defined.}
 \item{method}{the method for defining the threshold cycles. Either 'var' for a user-defined value or 'cpD2' (overrides 'var') for estimation from the second derivative maximum of the first run.}
 \item{dil}{a vector with the concentration (or dilution) values corresponding to the qPCR runs.}
 \item{pred}{a vector of threshold cycles whose concentrations are to be predicted.}
 \item{fct}{the model used for the threshold cycle estimation. Any sigmoidal model in the 'drc' package. Defaults to 'l5()'.} 
 \item{plot}{logical. Defaults to \code{TRUE}. If \code{FALSE}, the values are returned.}
 \item{conf.lev}{the p-value for the confidence interval. Defaults to 0.05, can be omitted with \code{NULL}.}
 \item{opt.min}{numeric. The minimal value for the regression optimization. Must be tweaked by the user.}
 \item{opt.plot}{logical. If \code{TRUE}, the Fluorescence vs. AIC plot is displayed (with the best value highlighted).}
}

\value{
  A list with the following components:
  \item{ct}{the calculated threshold cycles.}
  \item{r.squared}{The r-square of the linear fit.}
  \item{adj.r.squared}{The adjusted r-square of the linear fit.}
  \item{AIC}{the AIC value of the linear fit.}
  \item{Eff}{the qPCR efficiency as calculated from the slope of the dilution curve.}
  \item{pred}{the predicted values for the data given in 'pred'.}
  \item{thresh.opt}{the optimal threshold value.}
  \item{AIC.opt}{the minimal AIC of the regression optimization.}
}

\author{
  Andrej-Nikolai Spiess
}

\examples{
###define dilutions and values to be predicted###
dil <- c(100000, 10000, 1000, 100, 10, 1)
pred <- c(17.2, 21.2, 26.4)
###use reps data and single curves of each dilution###
###(i.e. F1.1, F2.1, F3.1 ...)
###threshold fluorescence is 0.2)
m1 <- calib(reps, f.col = c(2,6,10,14,18,22), groups = NULL, method = "var", 
	thresh = 0.2, dil = dil, pred = pred, plot = TRUE)
print(m1)
###use replicates, and threshold fluorescence calculated from the first replicate group### 
dil <- c(100000, 10000, 1000)
m2 <- calib(reps, f.col = 2:13, groups = c(1,1,1,1,2,2,2,2,3,3,3,3), method = "cpD2", 
	thresh = 0.2, dil = dil, pred = pred, plot = TRUE)
print(m2)
###optimize the threshold cycle###
\dontrun{
dil <- c(100000, 10000, 1000, 100, 10, 1)
m1 <- calib(reps, f.col = c(2,6,10,14,18,22), groups = NULL, method = "cpD2", 
        thresh = 0.2, dil = dil, plot = TRUE, opt.min = 0.05, opt.plot = TRUE)
}
}

\keyword{models}
\keyword{nonlinear}
