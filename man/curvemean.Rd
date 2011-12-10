\name{curvemean}
\alias{curvemean}

\title{Building a model which averages a batch of qPCR curves}

\description{
Starting with a batch of qPCR curves from class \code{modlist}, the function builds a sigmoidal model which averages
 the cycle numbers or expression values of the curves at every y-value (raw fluorescence). Thus, in contrast to existing qPCR
 averaging methods, not only the efficiencies or threshold cycles (as can be obtained from multiple reference
 genes) are averaged, but the complete structure of the batch at every fluorescence value.
 Beware: This is curve averaging, NOT model averaging!  
}

\usage{
curvemean(ml, mean = c("cmean1", "cmean2", "amean", "gmean", "hmean"), 
          which = 1, plot = TRUE)
}

\arguments{
  \item{ml}{a qPCR batch object of class 'modlist'.}
  \item{mean}{the averaging function, see 'Details'. Default is the expression mean.}     
  \item{which}{the position of the curve in the list which is used for defining the y-values. Should be the curve with the lowest threshold cycle.}  	
  \item{plot}{should results be plotted?}
}

\details{
Starting from the raw fluorescence values or expression values of a defined curve in the list, the cycle number of all other curves at this value are
 calculated and averaged with the method as under \code{mean}. After that, a new sigmoidal model is fit to the obtained data, using the same sigmoidal model as in \code{ml}.
 Note: often works better (more data points can be averaged) if \code{norm = TRUE} when buidling the model list with \code{\link{modlist}}. 
 Five different averaging functions can be used:

Expression mean with efficiency = 2 ("cmean1"):
\deqn{\bar x_{exp} = \frac{-log\left(\frac{\sum_{i=1}^n2^{-x_i}}{n}\right)}{log(2)}}

Expression mean with averaged efficiency from the runs ("cmean2"):
\deqn{\bar x_{exp} = \frac{-log\left(\frac{\sum_{i=1}^n\bar{E}^{-x_i}}{n}\right)}{log(\bar{E})}}

Arithmetic mean ("amean"):
\deqn{\bar x = \frac{\sum_{i=1}^nx_i}{n}} 

Geometric mean ("gmean"):
\deqn{\bar x = \left(\prod_{i=1}^nx_i \right)^\frac{1}{n}}

Harmonic mean ("hmean"):
\deqn{\bar x = \frac{n}{\sum_{i=1}^n\frac{1}{x_i}}}

Why is \code{mean = "cmean1"} the default? This is because all classical versions of the 'mean' do not work in the context of real-time PCR. Consider the following example: We have a qPCR experiment with three biological replicates, yielding threshold cycles 25, 26, and 27. This means (with an efficiency = 2) that the second sample has half (1/2^1) and the third sample 1/fourth (1/2^2) amount of initial copy numbers. If we plainly average (arithmetic mean) the threshold cycles, we obtain 26. But if we do this for the initial copy numbers, we obtain, for example, (100 + 50 + 25)/3 = 58.33. But ct = 26 corresponds to 50 copies, so the arithmetic mean of the threshold cycles significantly underestimates the initial copy numbers. The same accounts for the geometric and harmonic mean, as they all do not consider the exponential nature of the amplification. In contrast to this, the 'expression mean' does: cmean1 = -log(sum(2^-25 + 2^-26 + 2^-27)/3)/log(2) =  25.77761. And this values represents the averaged copy numbers because (calculated as ratios) 100/58.33 gives the same as 2^25.77761/2^25. Therefore, the expression mean calculates a cycle number which corresponds to the initial copy numbers of the replicates. If \code{mean = "cmean2"} is applied, the efficiency is calculated as the arithmetic mean of the efficiencies from all curves.
}


\value{
A new model of class 'nls' and 'pcrfit' from the averaged curve.  
}

\author{
Andrej-Nikolai Spiess
}

\references{
Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes.\cr
Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A & Speleman F.\cr
\emph{Genome Biology} (2002), \bold{3}: research0034.1-0034.11.\cr
}

\examples{
## create arithmetic mean curve of four
## serial dilutions
ml1 <- modlist(reps, fluo = c(2, 6, 10, 14))
curvemean(ml1, mean = "amean")

## effect of normalizing data to [0, 1]:
## more averaged datapoints
ml2 <- modlist(reps, fluo = c(2, 6, 10, 14), norm = TRUE)
curvemean(ml2, mean = "amean")

## averaging the expression value
## (i.e. copy numbers)
## using efficiency = 2
curvemean(ml2, mean = "cmean1")

## using averaged efficiencies from
## all curves at each cycle
curvemean(ml2, mean = "cmean2")
}

\keyword{models}
\keyword{nonlinear}
