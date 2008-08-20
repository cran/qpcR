\name{PRESS}
\alias{PRESS}

\title{Allen's PRESS (Prediction Sum-Of-Squares) statistic}

\description{
 Calculates the PRESS statistic, a leave-one-out refitting and prediction method, as described in Allen (1971).
 Works for any regression model with a \code{call} slot and a \code{\link{predict}} function, hence all models of class \code{lm}, \code{glm}, \code{nls} and \code{drc} (and maybe more...). 
}

\usage{
PRESS(model)
}

\arguments{
\item{model}{a fitted model of the classes as described above.}
}

\details{
From a fitted model, each of the predictors \eqn{x_i, i = 1 \ldots{n}} is removed and the model is refitted to the \eqn{n-1} points.
The predicted value \eqn{\hat{y}_{i, -i}} is calculated at the excluded point \eqn{x_i} and the PRESS statistic is given by:
\deqn{\sum_{i=1}^n (y_i - \hat{y}_{i, -i})^2}
The PRESS statistic is a surrogate measure of crossvalidation of small sample sizes and a measure for internal validity.
Small values indicate that the model is not overly sensitive to any single data point.
}

\value{
A list with the following components:
\item{stat}{The PRESS statistic.}
\item{resid}{a vector containing the PRESS residuals for each \eqn{x_i}.}  
}

\author{
Andrej-Nikolai Spiess
}

\note{
There is also a \code{PRESS} function in library 'MPV' that works solely for \code{lm} models using the hat matrix.
}

\references{
The relationship between variable selection and data augmentation and a method for prediction.\cr
Allen DM.\cr
Technometrics (1974), \bold{16}:25-127.

The Prediction Sum of Squares as a Criterion for Selecting Predictor Variables.\cr
Allen DM.\cr
Technical Report Number 23 (1971), Department of Statistics, University of Kentucky.\cr

Classical and Modern Regression with Applications.\cr
Myers RH.\cr
Second Edition (1990), Duxbury Press (PWS-KENT Publishing Company), 299-304.
}

\examples{
## example for PCR analysis
mod <- pcrfit(reps, 1, 2, l5())
PRESS(mod)$stat
## or
mod <- multdrc(F1.1 ~ Cycles, data = reps, fct = l5())
PRESS(mod)$stat

## compare PRESS statistic in models
## with fewer parameters
mod2 <- update(mod, fct = l4())
PRESS(mod2)$stat
mod3 <- update(mod, fct = l3())
PRESS(mod3)$stat

## example for linear regression
x <- 1:10
y <- rnorm(10, x, 0.1)
mod <- lm(y ~ x)
PRESS(mod)$stat

## example for NLS fitting
## starting values need to be attached!
DNase1 <- subset(DNase, Run == 1)
inits <- getInitial(density ~ SSlogis(log(conc), Asym, xmid, scal), data = DNase1)
attach(as.list(inits))
fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
res <- PRESS(fm1DNase1)

## PRESS residuals plot
plot(fitted(fm1DNase1), res$resid)















}

\keyword{models}
\keyword{nonlinear}
