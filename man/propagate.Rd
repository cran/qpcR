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
The data of the original dataset is shuffled \code{nperm} times (with replacement) according to the ties defined in \code{group}.
In detail, two datasets are created for each permutation:
Dataset1 samples with replacement the data in the columns, such that the data is not shuffled between columns. Dataset2 is obtained by
  sampling with replacement between the columns as defined in \code{group}. For both datasets, the rows are evaluated and statistics are collected.
The confidence interval is calculated from all evaluations of Dataset1. A p-value is calculated from all events mean(Dataset2) smaller/bigger mean(Dataset2) divided by \code{nperm}.
Thus, the p-value gives a measure against the null hypothesis that the result in the initial group is just by chance.
Omitting a grouping definition will result in shuffling the data between all columns (which rarely makes sense). If columns of
 the data are to be fixed in Dataset2, these have to be defined by \code{NA}'s. See 'Details' and 'Examples'.
 
3) \bold{Error propagation:}\cr
For all variables in the original dataset, mean and standard deviation are calculated. Through gaussian error propagation (first-order Taylor expansion),
 the propagated error is calculated using a matrix approach (See 'Details') by either omitting or including the dataset covariance structure.   
}

\usage{
propagate(expr, data, type = c("raw", "stat"), do.sim = TRUE, use.cov = FALSE, 
          nsim = 10000, do.perm = FALSE, perm.crit = "m1 < m2", ties = NULL, 
          nperm = 2000, alpha = 0.05, plot = TRUE, xlim = NULL, ...)  
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
  \item{perm.crit}{a character string defining the differences in means for the permutation p-value. See 'Details'.}
  \item{ties}{a vector defining the columns in which permutation and data swapping should be applied in order to analyze if differences in groups are only by chance. See 'Details'.}
  \item{nperm}{the number of permutations to be performed.}  
  \item{alpha}{the confidence level.}
  \item{plot}{logical. Should histograms with confidence intervals (in blue) be plotted for all analyses?}
  \item{xlim}{the x-axis limits for tweaking the histograms (in case of severe skewness).}
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
A high tendency from deviation of normality is encountered in formulae where the error of the denominator is relatively high
 or in exponential models where the exponent has a high error. This is one of the problems that is inherent in real-time PCR analysis, as the classical
 ratio calculation with efficiencies (i.e. by the delta-ct method) is usually of the exponential type. 
The \code{ties} define the swapping order to create the 'NULL' dataset from which permutation p-values are calculated against the initial dataset.
 Same numbers define the columns to swap against, as \code{NA} defined columns are only permutated within and not swapped.

Examples:
Randomly allocate values between first/third and second/fourth column: c(1, 2, 1, 2)\cr
Randomly allocate values between first/second column, the values of the remaining columns are not reallocated: c(1, 1, NA, NA).

The criterium for the p-value (\code{perm.crit}) has to be defined by the user.
For example, let's say we calculate some value 0.2 which is a ratio between two groups.
We would hypothesize that by randomly reallocating the values between the groups the mean values are not equal or smaller than in the initial data. We would thus define \code{perm.crit} as "m1 < m2" meaning that we want to test if the 
 mean of the initial data (m1) is frequently smaller than by the randomly allocated data (m2). 
}

\value{
  A list with the following components:       
  \item{mean.Sim}{the mean of all single row-wise evaluations from the Monte Carlo simulated data.}
  \item{sd.Sim}{the standard deviation of all single row-wise evaluations from the Monte Carlo simulated data.}
  \item{med.Sim}{the median of all single row-wise evaluations from the Monte Carlo simulated data.}
  \item{mad.Sim}{the median average deviation of all single row-wise evaluations from the Monte Carlo simulated data.}   
  \item{data.Sim}{the Monte Carlo simulated data.}
  \item{eval.Sim}{all single row-wise evaluations from the Monte Carlo simulated data.}
  \item{conf.Sim}{the confidence values for \code{eval.Sim}.}
  \item{mean.Perm}{the mean from all row-wise evaluations after dataset resampling (permutation).}
  \item{sd.Perm}{the standard deviation from all row-wise evaluations after dataset resampling (permutation).}
  \item{med.Perm}{the median from all row-wise evaluations after dataset resampling (permutation).}
  \item{mad.Perm}{the median average deviation from all row-wise evaluations after dataset resampling (permutation).}   
  \item{conf.Perm}{the confidence values obtained for all resamples.}
  \item{pval.Perm}{the permutation p-value between means of initial group permutations and group swapping permutations.}
  \item{eval.Prop}{the result from evaluating the expression with the mean values.}
  \item{error.Prop}{the propagated error.}
  \item{deriv.Prop}{a list containing the partial derivatives expressions for each variable.}
  \item{conf.Prop}{the confidence values of the propagated error, assuming normality (be cautious about this point!).}
  \item{covMat}{the covariance matrix used for the Monte-Carlo simulation and error propagation.}    
}

\author{
Andrej-Nikolai Spiess
}

\note{
A more elaborate description of the different error types can be found under \url{www.dr-spiess.de/qpcR/errors.pdf}.
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

\seealso{
Function \code{\link{ratiocalc}} for error analysis within qPCR ratio calculation.
}


\examples{
## from summary data just calculate 
## Monte-Carlo and propagated error.
EXPR <- expression(x/y)
x <- c(5, 0.1)
y <- c(1, 0.01)
DF <- cbind(x, y)
res <- propagate(expr = EXPR, data = DF, type = "stat")

## Do Shapiro-Wilks test on Monte Carlo errors 
## !maximum 5000 datapoints can be used!
## => p.value indicates normality
shapiro.test(res$eval.Sim[1:5000])
## How about a graphical analysis:
qqnorm(res$eval.Sim)

## using raw data
## bring all vectors to same 
## length with NA's
EXPR <- expression(x*y)
x <- c(2, 2.1, 2.2, 2, 2.3, 2.1)
y <- c(4, 4, 3.8, 4.1, NA, NA)
DF <- cbind(x, y)  
res <- propagate(EXPR, DF, type = "raw")

## For replicate data, using relative 
## quantification ratio from qPCR.
## How good is the estimation of the propagated error?
## Done without using covariance in the 
## calculation and simulation.
## STRONG deviation from normality!
EXPR <- expression((E1^cp1)/(E2^cp2))
E1 <- c(1.73, 1.75, 1.77)
cp1 <- c(25.77,26.14,26.33)
E2 <-  c(1.72,1.68,1.65)
cp2 <- c(33.84,34.04,33.33)
DF <- cbind(E1, cp1, E2, cp2)
res <- propagate(EXPR, DF, type = "raw", do.sim = TRUE)
par(ask = TRUE)
shapiro.test(res$eval.Sim[1:5000])
qqnorm(res$eval.Sim)

## same setup as above but also
## using a permutation approach
## for resampling confidence interval
## Similar to what REST software will do!
## Scenario 1: just swapping the 
## threshold cycles
\dontrun{
res <- propagate(EXPR, DF, type = "raw", do.sim = TRUE, 
                  do.perm = TRUE, ties = c(NA, 1, NA, 1)) 

                  
## 95\% confidence interval for the ratio
res$conf.Perm
}

## Using covariances in calculation and simulation
res2 <- propagate(EXPR, DF, type = "raw", 
                use.cov = TRUE, do.sim = TRUE)
              
## Proof that covariance of Monte-Carlo
## simulated dataset is the same as from 
## initial data
res2$covMat
cov(res2$data.Sim)

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
propagate(expression(m*7 + b), DF, type = "stat", use.cov = CM)

## In a x/y regime, when does the propagated error start to
## deviate from normality if error of denominator increases?
## Watch increasing skewness of histogram!
\dontrun{
x <- c(5, 1)
for (i in seq(0.01, 0.5, by = 0.01)) {
  y <- c(1, i)
  DF <- cbind(x, y)
  res  <-  propagate(expression(x/y), DF, type = "stat", 
                      do.sim = TRUE, plot = FALSE)
  hist(res$data.Sim, nclass = 100, main = paste("sd(y):",i))      
}
}

}

\keyword{distribution}
\keyword{htest}
