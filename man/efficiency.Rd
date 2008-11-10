\name{efficiency}
\alias{efficiency}

\title{Calculation of qPCR efficiency and several other important qPCR parameters}

\description{
  This function calculates the PCR efficiency of a fitted 'drc' or 'mchoice' model and additionally 
  several other important values for qPCR quantification like the first and second derivatives and the 
  corresponding maxima thereof (i.e. threshold cycles). These values can be subsequently used for the 
  calculation of PCR kinetics, fold induction etc. All values are included in a graphical output of the fit.
  Additionally, several measures of goodness-of-fit are calculated, i.e. the Akaike Information Criterion (AIC),
  the residual variance and the R-square value.
  For relative and absolute quantification, the initial fluorescence F(0) is calculated and can be used for
  optical calibration if the template amount is defined (see 'Details'). 
}

\usage{
  efficiency(object, plot = TRUE, type = "cpD2", shift = 0, amount = NULL)
}

\arguments{
  \item{object}{an object of class 'drc' or 'mchoice' model.}
  \item{plot}{logical. If TRUE, a graph is displayed. If FALSE, values are printed out.}
  \item{type}{the method of efficiency estimation. See 'Details'.}
  \item{shift}{a user defined shift in cycles from the values defined by \code{type}. See Examples.}
  \item{amount}{the template amount or molecule number for quantitative calibration.}
}

\details{
The efficiency is always calculated from the efficiency curve (in blue), which is calculated according to \eqn{E = \frac{F(n)}{F(n-1)}} from the 
fitted curve, but taken from different points at the curve: 

"cpD2"  taken from the maximum of the second derivative curve,\cr
"cpD1"  taken from the maximum of the first derivative curve,\cr
"maxE"  taken from the maximum of the efficiency curve,\cr
"expR"  taken from the exponential region by \eqn{expR = cpD2-(cpD1-cpD2))},\cr
"Cy0"   as calculated according to Guescini et al. or\cr
a numeric value (taken from the threshold cycle output of the PCR software, i.e. 15.24).

The initial fluorescence F(0) for relative or absolute quantification is calculated by \eqn{F(0) = \frac{F(n)}{E(n)^{Cyc}}}, 
 with F(n) = raw fluorescence at the cycle number defined by \code{type}, E(n) = PCR efficiency at the cycle number defined by \code{type} and
 Cyc = the cycle number defined by \code{type}. If a template amount is defined, a conversion factor \eqn{cf = \frac{amount}{F(0)}} is given.
The different measures for goodness-of-fit give an overview for the validity of the efficiency estimation. 
First and second derivatives are calculated from the fitted function and the maxima of the derivatives curve and the efficiency curve are obtained. 
This function calls the first derivatives from the 'drc' library and the internal functions \code{deriv2.b} and \code{deriv2.l} for the calculation of second derivatives.
}

\value{
  A list with the following components:
  \item{eff}{the PCR efficiency.}
  \item{resVar}{the residual variance.}
  \item{AICc}{the bias-corrected Akaike Information Criterion.}
  \item{AIC}{the Akaike Information Criterion.}
  \item{Rsq}{the R-square value.}
  \item{cpD1}{the first derivative maximum (point of inflection in 'l4' or 'b4' models, 
  can be used for defining the threshold cycle).}
  \item{cpD2}{the second derivative maximum (turning point of cpD1, more often used for defining 
  the threshold cycle).}
  \item{cpE}{the PCR cycle with the highest efficiency, if "maxE" method is chosen.}
  \item{cpR}{the PCR cycle within the exponential region calculated as under 'Details', if "expR" method is chosen.}
  \item{Cy0}{the PCR threshold cycle 'Cy0' according to Guescini et al. See 'Details'.} 
  \item{fluo}{the raw fluorescence value at the point defined by \code{type}.}
  \item{init}{the initial template fluorescence, i.e. F(0), calculated as under 'Details'.}
  \item{cf}{the conversion factor between raw fluorescence and template amount, if the latter is defined.}
}

\author{
  Andrej-Nikolai Spiess & Christian Ritz
}

\references{
  Weihong Liu and David A. Saint (2002) Validation of a quantitative method for real time PCR kinetics, 
  \emph{BBRC}, \bold{294}, 347 - 353.
  A new real-time PCR method to overcome significant quantitative inaccuracy due to slight amplification inhibition.\cr
  Guescini M et al, \emph{BMC Bioinformatics}, 2008, \bold{9}: 326.  
}

\note{
  Three parameter models ('b3' or 'l3') do not work very well in calculating the PCR efficiency. 
  It is advisable not to take too many cycles of the plateau phase prior to fitting the model as this has a strong 
  effect on the validity of the efficiency estimates. 
}

\examples{
## Fitting initial model
m1 <-  multdrc(F1.1 ~ Cycles, data = reps, fct = l4())
efficiency(m1)
 
## of best model
m2 <- mchoice(m1)
efficiency(m2) #uses type="cpD2"

## using "maxE" method, with calculation of PCR efficiency
## one cycle 'upstream' of the cycle of maximum efficiency
efficiency(m2, type = "maxE", shift = 1)

##using the exponential region
efficiency(m2, type = "expR")

## using threshold value (15.32) from PCR software
efficiency(m2, type = 15.32)
 
## using the first 30 cycles and a template amount
## (optical calibration)
m3 <-  multdrc(F1.1 ~ Cycles, data = reps[1:30,], fct = l5())
efficiency(m3, amount = 1E3)

}

\keyword{models}
\keyword{nonlinear}
