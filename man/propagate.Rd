\name{propagate}
\alias{propagate}
\encoding{latin1}

\title{General error analysis function using different methods}

\description{
A general function for the calculation of errors. Can be used for qPCR data, but any data that should be subjected to error analysis will do. The different error types can be calculated for any given expression from either replicates or from statistical summary data (mean & standard deviation). 
}

\usage{
propagate(expr, data, type = c("raw", "stat"), do.sim = FALSE, 
          dist.sim = c("norm", "tnorm"), use.cov = FALSE, nsim = 10000, 
          do.perm = FALSE, perm.crit = NULL, ties = NULL, nperm = 2000, 
          alpha = 0.05, plot = TRUE, logx = FALSE, verbose = FALSE, ...)  
}

\arguments{
  \item{expr}{an expression, such as \code{expression(x/y)}.}
  \item{data}{a dataframe or matrix containing either a) the replicates in columns or b) the means in the first row and the standard deviations in the second row. The variable names must be defined in the column headers.}
  \item{type}{either \code{raw} if replicates are given, or \code{stat} if means and standard deviations are supplied.}
  \item{do.sim}{logical. Should Monte Carlo error analysis be applied?}
  \item{dist.sim}{\code{"norm"} will use a multivariate normal distribution, \code{"tnorm"} a multivariate truncated normal distribution. See 'Details'.}
  \item{use.cov}{logical or variance-covariance matrix with the same column descriptions and column order as \code{data}. If \code{TRUE} together with replicates, the covariances are calculated from these and used within the Monte-Carlo simulation and error propagation. If \code{type = "stat"}, a square variance-covariance matrix can be supplied in the right dimensions (n x n, n = number of variables). If \code{FALSE}, the Monte-Carlo simulation and error propagation use only the diagonal (variances).}
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
The implemented methods are:\cr\cr
1) \bold{Monte-Carlo simulation:}\cr
For each variable in the dataset, simulated data with \code{nsim} samples is generated from a multivariate normal distribution using the mean \eqn{\mu} and standard deviation \eqn{\sigma} of each variable. All data is coerced into a new dataset that has the same covariance structure as the initial dataset. Each row of the simulated dataset is evaluated and statistics are calculated from the \code{nsim} evaluations. In scenarios that are nonlinear in nature (i.e. exponential setups) the distribution of the result values can be extremely skewed, mainly due to the simulated values at the extreme end of the normal distribution. Setting \code{dist.sim = "tnorm"} will fit a multivariate normal distribution,  calculate the lower/upper 2.5\% quantile on each side for each input variable and use these as bounds for simulating from a multivariate truncated normal distribution. This will (in part) remove some of the skewness in the result distribution. 

2) \bold{Permutation approach:}\cr
The data of the original dataset is permutated \code{nperm} times by binding observations together according to \code{ties}. The \code{ties} bind observations that can be independent measurements from the same sample. In qPCR terms, this would be a real-time PCR for two different genes on the same sample. If \code{ties} is omitted, the observations are shuffled independently. In detail, two datasets are created for each permutation: Dataset1 samples the rows (replicates) of the data according to \code{ties}. Dataset2 is obtained by sampling the columns (samples), also binding columns as in \code{ties}. For both datasets, the permutations are evaluated and statistics are collected. The confidence interval is calculated from all evaluations of Dataset1. A p-value is calculated from all permutations that follow \code{perm.crit}, whereby \code{init} reflects the permutations of the initial data and \code{perm} the permutations of the randomly reallocated samples. Thus, the p-value gives a measure against the null hypothesis that the result in the initial group is just by chance. See also 'Examples'.\cr
The criterium for the permutation p-value (\code{perm.crit}) has to be defined by the user. For example, let's say we calculate some value 0.2 which is a ratio between two groups. We would hypothesize that by randomly reallocating the values between the groups the mean values are not equal or smaller than in the initial data. We would thus define \code{perm.crit} as "perm < init" meaning that we want to test if the mean of the initial data (\code{init}) is frequently smaller than by the randomly allocated data (\code{perm}). The default (\code{NULL}) is to test all three variants "perm > init", "perm == init" and "perm < init".
 
3) \bold{Error propagation:}\cr
The propagated error is calculated by gaussian error propagation (first-order Taylor expansion) using a matrix approach. Often omitted, but important in models where the variables are dependent (i.e. linear regression), is the second covariance term.
\deqn{\sigma_Y^2 = \sum_{i}\left(\frac{\partial f}{\partial X_i} \right)^2 \sigma_i^2 + \sum_{i \neq j}\sum_{j \neq i}\left(\frac{\partial f}{\partial X_i} \right)\left(\frac{\partial f}{\partial X_j} \right) \sigma_{ij}} 
\code{propagate} calculates the propagated error either with or without covariances by using the matrix representation \deqn{C_Y = F_XC_XF_X^T} with \eqn{C_Y} = the propagated error, \eqn{F_X} = the p x n matrix with the results from the partial derivatives, \eqn{C_X} = the p x p variance-covariance matrix and \eqn{F_X^T} = the n x p transposed matrix of \eqn{F_X}. Depending on the input formula, the error propagation may result in an error that is not normally distributed. The Monte Carlo simulation, starting with normal distributions of the variables, can clarify this. The plots obtained from this function will also clarify this potential caveat. A high tendency from deviation of normality is encountered in formulas where the error of the denominator is relatively high or in exponential models where the exponent has a high error. This is one of the problems that is inherent in real-time PCR analysis, as the classical ratio calculation with efficiencies (i.e. by the delta-ct method) is usually of the exponential type.  
}

\value{
A plot containing the histograms of the Monte-Carlo simulation, the permutation values and histogram of the error propagation. Additionally inserted in the plots are a boxplot, the median values in red and the confidence intervals as blue borders. 

A list with the following components (if \code{verbose = TRUE}):   
  \item{data.Sim}{the Monte Carlo simulated data with the evaluations in the last column.}       
  \item{data.Perm}{the data of the permutated observations and samples with the corresponding evaluations and the decision according to \code{perm.crit}.}      
  \item{data.Prop}{\code{nsim} values taken from a normal distribution with mean and s.d. calculated from the propagated error.}   
  \item{derivs}{a list containing the partial derivatives expressions for each variable.}   
  \item{covMat}{the covariance matrix used for the Monte-Carlo simulation and error propagation.}    
  \item{summary}{a summary of the collected statistics, given as a dataframe. These are: mean, s.d. median, mad, lower/upper confidence interval and permutation p-value(s).} 
}

\author{
Andrej-Nikolai Spiess
}   

\references{
\bold{Error propagation (in general):}\cr
An Introduction to error analysis.\cr
Taylor JR.\cr
University Science Books (1996), New York.\cr

A very good technical paper describing error propagation by matrix calculation can be found under \url{www.nada.kth.se/~kai-a/papers/arrasTR-9801-R3.pdf}.\cr

Evaluation of measurement data - Guide to the expression of uncertainty in measurement.\cr
JCGM 100:2008 (GUM 1995 with minor corrections).\cr
\url{http://www.ifcc.org/pdf/GUM_JCGM_100_2008_E.pdf}.\cr

\bold{Error propagation (in qPCR):}\cr
Error propagation in relative real-time reverse transcription polymerase chain reaction quantification models: The balance between accuracy and precision.\cr
Nordgard O, Kvaloy JT, Farmen RK, Heikkila R.\cr
\emph{Analytical Biochemistry} (2006), \bold{356}: 182-193.\cr

qBase relative quantification framework and software for management and analysis of real-time quantitative PCR data.\cr
Hellemans J, Mortier G, De Paepe A, Speleman F, Vandesompele J.\cr
\emph{Genome Biology} (2007), \bold{8}: R19.\cr 

\bold{Multivariate normal distribution:}\cr
Stochastic Simulation.\cr
Ripley BD.\cr
Stochastic Simulation (1987). Wiley. Page 98.\cr

\bold{Testing for normal distribution:}\cr
Testing for  Normality.\cr
Thode Jr. HC.\cr
Marcel Dekker (2002), New York.\cr

Approximating the Shapiro-Wilk W-test for non-normality.\cr
Royston P.\cr 
\emph{Statistics and Computing} (1992), \bold{2}: 117-119.\cr
}

\seealso{
Function \code{\link{ratiocalc}} for error analysis within qPCR ratio calculation.
}

\examples{
## From summary data just calculate 
## Monte-Carlo and propagated error.
EXPR <- expression(x/y)
x <- c(5, 0.01)
y <- c(1, 0.01)
DF <- cbind(x, y)
RES1 <- propagate(expr = EXPR, data = DF, type = "stat", 
                 do.sim = TRUE, verbose = TRUE)

## Do Shapiro-Wilks test on Monte Carlo evaluations 
## !maximum 5000 datapoints can be used!
## => p.value on border to non-normality
shapiro.test(RES1$data.Sim[1:5000, 3])
## How about a graphical analysis:
qqnorm(RES1$data.Sim[, 3])

## Using raw data
## If data is of unequal length,
## use qpcR:::cbind.na to avoid replication!
## Do permutations (swap x and y values)
## and simulations.
EXPR <- expression(x*y)
x <- c(2, 2.1, 2.2, 2, 2.3, 2.1)
y <- c(4, 4, 3.8, 4.1, 3.1)
DF <- qpcR:::cbind.na(x, y)  
RES2 <- propagate(EXPR, DF, type = "raw", do.perm = TRUE, 
                 do.sim = TRUE, verbose = TRUE)
RES2$summary

## Example in Annex H.1 from the GUM manual
## (see 'References'), an end gauge calibration
## study
EXPR <- expression(ls + d - ls * (da * the + as * dt))
ls <- c(50000623, 25)
d <- c(215, 9.7)
da <- c(0, 0.58E-6)
the <- c(-0.1, 0.41)
as <- c(11.5E-6, 1.2E-6)
dt <- c(0, 0.029)
DF <- cbind(ls, d, da, the, as, dt)
propagate(expr = EXPR, data = DF, type = "stat", 
                    do.sim = TRUE, verbose = TRUE)
## => u = 31.71 (GUM says 32).  

## For replicate data, using relative 
## quantification ratio from qPCR.
## How good is the estimation of the propagated error?
## Done without using covariance in the 
## calculation and simulation.
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
RES3 <- propagate(EXPR, DF, type = "raw", do.sim = TRUE,
                 do.perm = FALSE, verbose = TRUE, logx = TRUE)                 
## STRONG deviation from normality!
shapiro.test(RES3$data.Sim[1:5000, 5])
qqnorm(RES3$data.Sim[, 5])

## Same setup as above but also
## using a permutation approach
## for resampling the confidence interval.
## Cp's and efficiencies are tied together
## because they are two observations on the
## same sample! 
## Similar to what REST2008 software does.
RES4 <- propagate(EXPR, DF, type = "raw", do.sim = TRUE,
                 perm.crit = NULL, do.perm = TRUE, 
                 ties = c(1, 1, 2, 2), logx = TRUE, verbose = TRUE)
RES4$summary
## p-value of 0 in perm < init indicates that not a single 
## exchange of group memberships resulted in a smaller ratio!
              
## Proof that covariance of Monte-Carlo
## simulated dataset is the same as from 
## initial data
RES4$covMat
cov(RES4$data.Sim[, 1:4])
all.equal(RES4$covMat, cov(RES4$data.Sim[, 1:4]))

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
RES5 <- propagate(expression(m*7 + b), DF, type = "stat", use.cov = CM)
RES5

## In a x/y regime, when does the propagated error start to
## deviate from normality if error of denominator increases?
## Watch increasing deviation in qqnorm-plot!
\dontrun{
x <- c(5, 1)
for (i in seq(0, 0.2, by = 0.01)) {
  y <- c(1, i)
  DF <- cbind(x, y)
  RES6  <-  propagate(expression(x/y), DF, type = "stat", 
                      do.sim = TRUE, plot = FALSE, verbose = TRUE,
                      logx = TRUE)
  qqnorm(RES6$data.Sim[, 3])        
}
}
}   

\keyword{distribution}
\keyword{htest}
