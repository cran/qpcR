\name{propagate}
\alias{propagate}
\encoding{latin1}

\title{General error analysis function using Monte-Carlo simulation/permutation/(first-order) error propagation}

\description{
A general function for the calculation of errors. Can be used for qPCR data, but any data that should be subjected to error analysis will do.
The different error types can be calculated for any given expression from either replicates or from statistical summary data (mean & standard deviation). 
These are:\cr\cr
1) \bold{Monte-Carlo simulation:}\cr
For each variable in the dataset, simulated data with \code{nsim} samples is generated from a multivariate normal distribution using
 the mean and s.d. of each variable. All data is coerced into a new dataset that has the same covariance structure as the initial dataset.
Each row of the simulated dataset is evaluated and statistics are calculated from the \code{nsim} evaluations.

2) \bold{Permutation approach:}\cr
The data of the original dataset is permutated \code{nperm} times by binding observations together according to \code{ties}.
 The \code{ties} bind observations that can be independent measurements from the same sample. In qPCR terms, this would be a real-time PCR for two different genes on the same sample. If \code{ties} is omitted, the observations are shuffled independently. 
In detail, two datasets are created for each permutation:
Dataset1 samples the rows (replicates) of the data according to \code{ties}. Dataset2 is obtained by sampling the columns (samples), also binding columns as in \code{ties}. 
For both datasets, the permutations are evaluated and statistics are collected.
The confidence interval is calculated from all evaluations of Dataset1. A p-value is calculated from all permutations that follow \code{perm.crit}, 
 whereby \code{init} reflects the permutations of the initial data and \code{perm} the permutations of the randomly reallocated samples.
Thus, the p-value gives a measure against the null hypothesis that the result in the initial group is just by chance.
See 'Details' and 'Examples'.
 
3) \bold{Error propagation:}\cr
For all variables in the original dataset, mean and standard deviation are calculated. Through gaussian error propagation (first-order Taylor expansion),
 the propagated error is calculated using a matrix approach (see 'Details') by either omitting or including the dataset covariance structure.   
}

\usage{
propagate(expr, data, type = c("raw", "stat"), do.sim = FALSE, 
          use.cov = FALSE, nsim = 10000, do.perm = FALSE, 
          perm.crit = "perm > init", ties = NULL, nperm = 2000, 
          alpha = 0.05, plot = TRUE, logx = FALSE, verbose = FALSE, ...)  
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{data}{a dataframe or matrix containing either a) the replicates in columns or b) the means in the first row and the standard deviations
	      in the second row. The variable names must be defined in the column headers.}
  \item{type}{either \code{raw} if replicates are given, or \code{stat} if means and standard deviations are supplied.}
  \item{do.sim}{logical. Should Monte Carlo error analysis be applied?}
  \item{use.cov}{logical or variance-covariance matrix with the same column descriptions and column order as \code{data}. If \code{TRUE} together with replicates, 
	     the covariances are calculated from these and used within the Monte-Carlo simulation and error propagation. If \code{type = "stat"}, a square variance-covariance matrix can be supplied in the right dimensions 
             (n x n, n = number of variables). If \code{FALSE}, the Monte-Carlo simulation and error propagation use only the diagonal (variances).}
  \item{nsim}{the number of simulations to be performed, minimum is 5000.}  
  \item{do.perm}{logical. Should permutation error analysis be applied?}     
  \item{perm.crit}{a character string of one or more criteria defining the null hypothesis for the permutation p-value. See 'Details'.}
  \item{ties}{a vector defining the columns that should be tied together for the permutations. See 'Details'.}
  \item{nperm}{the number of permutations to be performed.}  
  \item{alpha}{the confidence level.}
  \item{plot}{logical. Should histograms with confidence intervals (in blue) be plotted for all analyses?}
  \item{logx}{logical. Should the x-axis of the graphs have logarithmic scale?}
  \item{verbose}{logical. If \code{TRUE}, a longer output is given including the simulated data, derivatives, covariance matrix etc.}
  \item{...}{other parameters to be supplied to \code{hist}, \code{boxplot} or \code{abline}.}
}

\details{
The propagated error is calculated by gaussian error propagation. Often omitted, but important in models where the variables are dependent (i.e. linear regression),
 is the second covariance term.
\deqn{\sigma_Y^2 = \sum_{i}\left(\frac{\partial f}{\partial X_i} \right)^2 \sigma_i^2 + \sum_{i \neq j}\sum_{j \neq i}\left(\frac{\partial f}{\partial X_i} \right)\left(\frac{\partial f}{\partial X_j} \right) \sigma_{ij}} 
\code{propagate} calculates the propagated error either with or without covariances by using the matrix representation
\deqn{C_Y = F_XC_XF_X^T}
with \eqn{C_Y} = the propagated error, \eqn{F_X} = the p x n matrix with the results from the partial derivatives, \eqn{C_X} = the p x p variance-covariance matrix and
 \eqn{F_X^T} = the n x p transposed matrix of \eqn{F_X}.
Depending on the input formula, the error propagation may result in an error that is not normally distributed. The Monte Carlo simulation, starting with normal distributions
 of the variables, can clarify this. The plots obtained from this function will also clarify this potential caveat.  
A high tendency from deviation of normality is encountered in formulas where the error of the denominator is relatively high
 or in exponential models where the exponent has a high error. This is one of the problems that is inherent in real-time PCR analysis, as the classical
 ratio calculation with efficiencies (i.e. by the delta-ct method) is usually of the exponential type. 

The criterium for the permutation p-value (\code{perm.crit}) has to be defined by the user.
For example, let's say we calculate some value 0.2 which is a ratio between two groups.
We would hypothesize that by randomly reallocating the values between the groups the mean values are not equal or smaller than in the initial data. 
We would thus define \code{perm.crit} as "perm < init" meaning that we want to test if the mean of the initial data (\code{init}) is frequently smaller than by the randomly allocated data (\code{perm}).
One can also supply something like \code{perm.crit = c("perm > init", "perm == init", "perm < init")} to test for all combinations. 
}

\value{
A plot containing the histograms of the Monte-Carlo simulation, the permutation values and histogram of the error propagation. 
Additionally inserted in the plots are a boxplot, the median values in red and the confidence intervals as blue borders. 

A list with the following components:       
  \item{summary}{a summary of the collected statistics, given as a dataframe. These are: mean, s.d. median, mad, lower/upper confidence interval and permutation p-value(s).}  
  \item{data.Sim}{the Monte Carlo simulated data with the evaluations in the last column, if \code{verbose = TRUE}.}       
  \item{data.Perm}{the data of the permutated observations and samples with the corresponding evaluations and the decision according to \code{perm.crit}, if \code{verbose = TRUE}. See 'Examples'.}      
  \item{derivs}{a list containing the partial derivatives expressions for each variable, if \code{verbose = TRUE}.}   
  \item{covMat}{the covariance matrix used for the Monte-Carlo simulation and error propagation, if \code{verbose = TRUE}.}    
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
Royston P (1992). Approximating the Shapiro-Wilk W-test for non-normality.\cr
\emph{Statistics and Computing}, \bold{2}, 117-119.\cr
}

\seealso{
Function \code{\link{ratiocalc}} for error analysis within qPCR ratio calculation.
}

\examples{
## From summary data just calculate 
## Monte-Carlo and propagated error.
EXPR <- expression(x/y)
x <- c(5, 0.1)
y <- c(1, 0.01)
DF <- cbind(x, y)
res <- propagate(expr = EXPR, data = DF, type = "stat", 
                 do.sim = TRUE, verbose = TRUE)

## Do Shapiro-Wilks test on Monte Carlo evaluations 
## !maximum 5000 datapoints can be used!
## => p.value indicates normality
shapiro.test(res$data.Sim[, 3][1:5000])
## How about a graphical analysis:
qqnorm(res$data.Sim[, 3])

## Using raw data
## bring all vectors to same 
## length with NA's.
## Do permutations (swap x and y values)
## and simulations
EXPR <- expression(x*y)
x <- c(2, 2.1, 2.2, 2, 2.3, 2.1)
y <- c(4, 4, 3.8, 4.1, 3.1, 3)
DF <- cbind(x, y)  
res2 <- propagate(EXPR, DF, type = "raw", do.perm = TRUE, 
                 do.sim = TRUE, verbose = TRUE)
## Have a look at the results
res2$summary

## For replicate data, using relative 
## quantification ratio from qPCR.
## How good is the estimation of the propagated error?
## Done without using covariance in the 
## calculation and simulation.
## STRONG deviation from normality!
## cp's and efficiencies are tied together
## because they are two observations on the
## same sample!
## As we are using an exponential type function,
## better to logarithmize the x-axis.
EXPR <- expression((E1^cp1)/(E2^cp2))
E1 <- c(1.73, 1.75, 1.77)
cp1 <- c(25.77,26.14,26.33)
E2 <-  c(1.72,1.68,1.65)
cp2 <- c(33.84,34.04,33.33)
DF <- cbind(E1, cp1, E2, cp2)
res3 <- propagate(EXPR, DF, type = "raw", do.sim = TRUE,
                 do.perm = FALSE, verbose = TRUE, logx = TRUE)                 
par(ask = TRUE)
shapiro.test(res3$data.Sim[, 5][1:5000])
qqnorm(res3$data.Sim[, 5])

## Same setup as above but also
## using a permutation approach
## for resampling confidence interval
## Cp's and efficiencies are tied together
## because they are two observations on the
## same sample! 
## Similar to what REST2008 software does.
res4 <- propagate(EXPR, DF, type = "raw", do.sim = TRUE,
                 perm.crit = "perm < init", do.perm = TRUE, 
                 ties = c(1, 1, 2, 2), logx = TRUE, verbose = TRUE)
res4$summary
              
## Proof that covariance of Monte-Carlo
## simulated dataset is the same as from 
## initial data
res4$covMat
cov(res4$data.Sim[, 1:4])

## Error propagation in a linear model 
## using the covariance matrix from summary.lm
## Estimate error of y for x = 7
x <- 1:10	
set.seed(123)
y <- x + rnorm(10, 0, 1) ##generate random data	
mod <- lm(y ~ x)
summ <- summary(mod)
## make matrix of parameter estimates and standard error
DF <- t(coef(summ)[, 1:2]) 
colnames(DF) <- c("b", "m")
CM <- vcov(mod) ## take var-cov matrix
colnames(CM) <- c("b", "m")
res5 <- propagate(expression(m*7 + b), DF, type = "stat", use.cov = CM)
res5

## In a x/y regime, when does the propagated error start to
## deviate from normality if error of denominator increases?
## Watch increasing skewness of histogram!
\dontrun{
x <- c(5, 1)
for (i in seq(0.01, 0.5, by = 0.01)) {
  y <- c(1, i)
  DF <- cbind(x, y)
  res  <-  propagate(expression(x/y), DF, type = "stat", 
                      do.sim = TRUE, plot = FALSE, verbose = TRUE)
  hist(res$data.Sim[, 3], nclass = 100, main = paste("sd(y):",i))      
}
}
}   

\keyword{distribution}
\keyword{htest}
