\name{curvemean}
\alias{curvemean}

\title{Building a model which averages a batch of qPCR curves}

\description{
Starting with a batch of qPCR curves from class \code{modlist}, the function builds a sigmoidal model which averages
 the cycle numbers of the curves at every y-value (i.e. raw fluorescence). Thus, in contrast to existing qPCR
 averaging methods, not only the efficiencies or threshold cycles (as can be obtained from multiple reference
 genes) are averaged, but the complete structure of the batch at every fluorescence value.
 Beware: This is curve averaging, NOT model averaging!  
}

\usage{
  curvemean(modlist, mean = c("amean", "gmean", "hmean"), 
	    fct = l5(), which = 1)
}

\arguments{
  \item{modlist}{a qPCR batch object obtained from \code{\link{modlist}}.}
  \item{mean}{the averaging function, see 'Details'. Default is the arithmetic mean.}
  \item{fct}{the sigmoidal model to be used for the averaged curve.}  
  \item{which}{The position of the curve in the list which is used for defining the y-values. Defaults to the first curve.}  	
}

\details{
Starting from the raw fluorescence values of a defined curve in the list, the cycle number of all other curves at this value are
 calculated and averaged. After that, a new sigmoidal model is fit to the obtained data. 
 Note: often works better (more data points can be averaged) if \code{norm = TRUE} when buidling the model list with \code{\link{modlist}}. 
 Three different averaging functions can be used:

Arithmetic mean:
\deqn{\bar x = \frac{\sum_{i=1}^nx_i}{n}} 

Geometric mean:
\deqn{\bar x = \left(\prod_{i=1}^nx_i \right)^\frac{1}{n}}

Harmonic mean:
\deqn{\bar x = \frac{n}{\sum_{i=1}^n\frac{1}{x_i}}}
}

\value{
A model of class \code{drc} from the averaged curve.  
}

\author{
Andrej-Nikolai Spiess
}

\references{
Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes.
Vandesompele et al., \emph{Genome Biology}, 2002, \bold{3}: research0034.1-0034.11.
}

\examples{
## create arithmetic mean curve of four
## serial dilutions
ml <- modlist(reps, 1, c(2, 6, 10, 14))
curvemean(ml)

## effect of normalizing data to [0, 1]:
## more averaged datapoints
ml <- modlist(reps, 1, c(2, 6, 10, 14), norm = TRUE)
curvemean(ml)

## difference with geometric mean
curvemean(ml, mean = "gmean")
}

\keyword{models}
\keyword{nonlinear}
