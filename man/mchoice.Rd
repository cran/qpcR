\name{mchoice}
\alias{mchoice}

\encoding{latin1}

\title{Selection of the best model by nested F-tests/likelihood ratios/Akaike weights}

\description{
  Model selection by comparison of different models using\cr 

  1) the maximum log likelihood value,\cr
  2) Akaike's Information Criterion (AIC),\cr
  3) bias-corrected Akaike's Information Criterion (AICc),\cr
  4) the estimated residual variance,\cr
  5) the p-value from a nested F-test on the residual variance,\cr
  6) the p-value from the likelihood ratio (chi-square),\cr
  7) the Akaike weights based on AIC and\cr
  8) the Akaike weights based on AICc.

  The best model is chosen by 5), 6) or 7) and returned as a new model.
}

\usage{
  mchoice(object, fctList = NULL, sig.level = 0.05, verbose = TRUE, 
	  crit = c("ftest", "ratio", "weights"))
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{fctList}{a list of functions to be analyzed, i.e. for a non-nested regime. Should also contain the original model.}
  \item{sig.level}{the significance level for the nested F-test.}
  \item{verbose}{logical. If \code{TRUE}, the result matrix is displayed in the console.}
  \item{crit}{the criterium for model selection. Either 'ftest'/'ratio' for nested models or 'weights' for nested and non-nested models.}
 }

\details{
  Criteria 5) and 6) cannot be used for comparison unless the models are nested. 
  Criterion 7), Akaike weights, can be used for nested and non-nested regimes.
  For criterion 1) the larger the better. For criteria 2), 3) and 4): the smaller the better.
  The best model is chosen either from the nested F-test, likelihood ratio or Akaike weights and returned as a new model.
  When using 'ftest'/'ratio' the corresponding nested function are analyzed automatically, i.e. b3/b4/b5; l3/l4/l5; w3/w4.  
}

\value{
  A model of the best fit selected by the nested F-tests, likelihood ratios or Akaike weights.
  The new model has an additional list item 'retMat' with the result matrix from the criterion tests.
}

\author{
  Andrej-Nikolai Spiess
}

\seealso{
 \code{\link{LR}}, \code{\link{akaike.weights}}
}


\examples{
m1 <- multdrc(F1.1 ~ Cycles, data = reps, fct = l3())
### choose best model based on F-tests 
### on the corresponding nested models
m2 <- mchoice(m1)
summary(m2)  ### Converted to l5() model !
plot(m2, log = "")

### use Akaike weights on non-nested models
### compare to original model
m2 <- mchoice(m1, fctList = list(l3(), l5(), b3(), w4(), baro5()), crit = "weights")
}

\keyword{models}
\keyword{nonlinear}
