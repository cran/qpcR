\name{mchoice}
\alias{mchoice}

\title{Selection of the best model by nested F-tests}

\description{
  Model selection by comparison of different models using\cr 

  1) the maximum log likelihood value,\cr
  2) Akaike's Information Criterion,\cr
  3) bias-corrected Akaike's Information Criterion,\cr
  4) the estimated residual variance and\cr
  5) the p-value from a nested F-test on the residual variance.\cr

  The best model is chosen by 5) and returned as a new model.
}

\usage{
  mchoice(object, sig.level = 0.05, verbose = TRUE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{sig.level}{the significance level for the nested F-test.}
  \item{verbose}{logical. If \code{TRUE}, the result matrix is displayed in the console.}
 }

\details{
  Criterion 5) cannot be used for comparison unless the models are nested. 
  For criterion 1) the larger the better. For criteria 2), 3) and 4): the smaller the better.
  The best model is chosen from the nested F-test and returned as a new 'drc' model.
}

\value{
  A 'drc' model of the best fit selected by the nested F-tests.
  The new model has an additional list item 'retMat' with the result matrix from the criterion tests.
}

\author{
  Andrej-Nikolai Spiess & Christian Ritz
}


\examples{
m1 <- multdrc(F1.1 ~ Cycles, data = reps, fct = l3())
### choose best model
m2 <- mchoice(m1)
summary(m2)  ### Converted to l5() model !
plot(m2, log = "")
}

\keyword{models}
\keyword{nonlinear}
