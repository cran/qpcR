\name{propagate}
\alias{propagate}

\title{General error propagation function using covariates and Monte Carlo simulations}

\description{
A general function for calculating the propagated error. This can be calculated for any given function from either replicates or from statistical
 summary data (mean & standard deviation). If independency of variables cannot be assumed, covariates can be calculated from the replicated data
 or by supplying a variance-covariance matrix to the input. The resulting distribution of the propagated error can be evaluated by Monte Carlo simulation
 of the input by generating at least 5000 datapoints from the (normal or multivariate normal) distribution of the variables.     
}

\usage{propagate(fun, vals, type = c("raw", "stat"), cov = FALSE, do.sim = FALSE, 
		cov.sim = FALSE, nsim = 10000)  
}

\arguments{
  \item{fun}{a raw function description such as \code{x/y}, without any quotes.}
  \item{vals}{a dataframe or matrix containing either a) the replicates in columns or b) the means in the first row and the standard deviations
	      in the second row. The variable names must be defined in the column headers.}
  \item{type}{Either \code{raw} if replicates are given, or \code{stat} if means and standard deviations are supplied.}
  \item{cov}{logical or variance-covariance matrix with the same column descriptions and column order as \code{vals}. If \code{TRUE} together with replicates, 
	     the covariances are calculated from these. If \code{type = "stat"}, a square variance-covariance matrix can be supplied in the right dimensions 
             (n x n, n = number of variables). If \code{FALSE}, the error is propagated using only the variances.}
  \item{do.sim}{logical. Should Monte Carlo simulation be applied to the propagated error?}
  \item{cov.sim}{logical. If \code{TRUE}, the simulated data is generated with a multivariate covariance structure similar to the input variables.}
  \item{nsim}{the number of simulations to be performed, minimum is 5000.}  
}

\details{
The propagated error is calculated by gaussian error propagation. Often omitted, but important in models where the variables are dependent (i.e. linear regression),
 is the second covariance term.
\deqn{\sigma_Y^2 = \sum_{i}\left(\frac{\partial f}{\partial X_i} \right)^2 \sigma_i^2 + \sum_{i \neq j}\sum_{j \neq i}\left(\frac{\partial f}{\partial X_i} \right)\left(\frac{\partial f}{\partial X_j} \right) \sigma_{ij}} 
\code{propagate} calculates the propagated error either with or without covariances by using the matrix representation
\deqn{C_Y = F_XC_XF_X^T}
with \eqn{C_Y} = the propagated error, \eqn{F_X} = the p x n matrix with the results from the partial derivatives, \eqn{C_X} = the p x p variance-covariance matrix and
 \eqn{F_X^T} = the n x p transposed matrix of \eqn{F_X}.
Depending on the input formula, the error propagation may result in an error that is not normally distributed. The Monte Carlo simulation, starting with the distributions
 of the variables, can clarify this. For acquisition of even more realistic estimates from the simulation process, a multivariate dataset having the same
 covariance structure can be created. A high tendency from deviation of normality is encountered in formulae where the error of the denominator is relatively high
 or in exponential models where the base has a high error. This is one of the problems that is inherent in real-time PCR analysis, as the classical
 ratio calculation with efficiencies (i.e. by the delta-ct method) is usually of the exponential type. The user of this function will find that within ratio calculations
 the propagated error is often non-normal (see 'Examples'). This package uses essentially the same algorithm as \code{delta.method} in package 'alr3' but
 has more error and compatibility checking and univariate/multivariate simulation.  
}

\value{
  A list with the following components:
  \item{errProp}{the propagated error, obtained from the input variables.}
  \item{evalProp}{the result from the evaluated function, obtained from the input variables.}
  \item{errSim}{the average of the propagated errors obtained from the Monte Carlo simulation.}  
  \item{evalSim}{the average of the evaluated function obtained from the Monte Carlo simulation.}  
  \item{derivs}{a list with each of the partial derivatives as items.}
  \item{covMat}{the square variance-covariance matrix that was used for the calculations.}  
  \item{simVec}{the propagated errors from the simulation coerced into a vector with length \code{nsim}. To be used for subsequent analysis, i.e tests for normality.}
}

\author{
  Andrej-Nikolai Spiess
}

\references{
Error propagation (in general):\cr
Taylor JR (1996). An Introduction to error analysis. University Science Books, New York.\cr 
A very good technical paper describing error propagation by matrix calculation can be found under \url{www.nada.kth.se/~kai-a/papers/arrasTR-9801-R3.pdf}.\cr

Error propagation (in qPCR):\cr
Nordgard O \emph{et al.} (2006). Error propagation in relative real-time reverse transcription polymerase chain reaction quantification models: The balance between accuracy and precision. \emph{Analytical Biochemistry}, \bold{356}, 182-193.\cr
Hellemans J \emph{et al.} (2007). qBase relative quantification framework and software for management and analysis of real-time quantitative PCR data. \emph{Genome Biology}, \bold{8}: R19.\cr 

Multivariate normal distribution:\cr
Ripley BD (1987). Stochastic Simulation. Wiley. Page 98.\cr

Testing for normal distribution:\cr
Thode Jr. HC (2002). Testing for  Normality. Marcel Dekker, New York.\cr
Royston P (1992). Approximating the Shapiro-Wilk W-test for non-normality. \emph{Statistics and Computing}, \bold{2}, 117-119.\cr
}


\examples{
## Starting from summary data just calculate the propagated error.
x <- c(5, 0.1)
y <- c(1, 0.01)
DF <- cbind(x, y)
propagate(x/y, vals = DF, type = "stat")

## Do simulations and evaluate error distribution.
res <- propagate(x/y, vals = DF, type = "stat", cov = TRUE, do.sim = TRUE)
## Do Shapiro-Wilks test on simulated errors 
## !maximum 5000 datapoints can be used!
## => p.value indicates normality
shapiro.test(res$simVec[1:5000])
## How about a graphical analysis:
qqnorm(res$simVec)
## histogram also indicates no skewness or kurtosis
hist(res$simVec, nclass = 100)

## For replicate data, using relative quantification ratio from qPCR.
## How good is the estimation of the propagated error?
## Done without using covariance in the calculation and simulation.
## STRONG deviation from normality!
E1 <- c(1.73, 1.75, 1.77)
cp1 <- c(25.77,26.14,26.33)
E2 <-  c(1.72,1.68,1.65)
cp2 <- c(33.84,34.04,33.33)
DF <- cbind(E1, cp1, E2, cp2)
res <- propagate((E1^cp1)/(E2^cp2), vals = DF, type = "raw", do.sim = TRUE)
shapiro.test(res$simVec[1:5000])
qqnorm(res$simVec)

## Using covariances in calculation and simulation
res <- propagate((E1^cp1)/(E2^cp2), vals = DF, type = "raw", cov = TRUE, do.sim = TRUE, 
		cov.sim = TRUE)

## Error propagation in a linear model using the covariance matrix from summary.lm
## Estimate error of y for x=7
x <- 1:10	
set.seed(123)
y <- x + rnorm(10, 0, 1) ##generate random data	
mod <- lm(y ~ x)
summ <- summary(mod)
DF <- t(coef(summ)[, 1:2]) ## make matrix of parameter estimates and standard error
colnames(DF) <- c("b", "m")
CM <- vcov(mod) ## take var-cov matrix
colnames(CM) <- c("b", "m")
propagate(m*7 + b, vals = DF, type = "stat", cov = CM)

## In a x/y regime, when does the propagated error start to
## deviate from normality if error of denominator increases?
## Watch skewness of histogram!
\dontrun{
x <- c(5, 1)
for (i in seq(0.001, 0.2, by = 0.001)) {
      y <- c(1, i)
      DF <- cbind(x, y)
      res  <-  propagate(x/y, vals = DF, type = "stat", do.sim = TRUE)
      hist(res$simVec, nclass = 100, main = paste("sd(y):",i))      
}
}

}

\keyword{distribution}
\keyword{htest}
