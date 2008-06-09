\name{ratioplot}
\alias{ratioplot}

\title{Plot method for data of class 'ratiocalc' which displays barplots of ratios with error bars}

\description{
A plot method in which data obtained from \code{\link{ratiocalc}} is displayed as a barplot or hanging barplot. Error bars can be added to the bars
 that had been calculated by any of the four different error values of \code{\link{propagate}}, i.e.
 propagated/Monte Carlo propagated/evaluated/Monte Carlo evaluated errors. Either significance codes or values from the t-tests can be added to
 the bars.
}

\usage{ratioplot(data, errbar = c("prop", "propSim", "eval", "evalSim", "none"), 
	  sem = FALSE, order = NA, y.fac = 1.5, normcol = NA, errwid = NULL, 
	  type = c("bar", "hbar"), plot.t = c("none", "stars", "values"), 
	  offset = 1, ...)
}

\arguments{
  \item{data}{data of class 'ratiocalc'.}
  \item{errbar}{type of error that is used in displaying the error bars. See 'Details'.}
  \item{sem}{logical. Should the standard error of the mean (s.d./sqrt(n)) be displayed?}
  \item{order}{a numeric vector defining subsets or reordering of the plotted bars.}
  \item{y.fac}{multiplication factor for the y-axis extension, can be tweaked for a more appealing output.}
  \item{normcol}{column number that all other data should be normalized against. See 'Details'.}
  \item{errwid}{The width of the error bars. If \code{NULL} (default), calculated from the width of the bars.}
  \item{type}{the plot type. See 'Details'.}
  \item{plot.t}{how to display the results from the t-tests. Either significance codes (i.e. < 0.01 => "**") or the p-values.}
  \item{offset}{distance of \code{plot.t} from the end of the bars. For tweaking.} 
  \item{...}{other parameters to be passed to \code{\link{barplot}} or \code{\link{arrows}}.}
}

\details{
Ratios are displayed from all permutations/combinations of sample/replicate sample PCR runs, as described in \code{\link{ratiocalc}}. The data can also be normalized
 against one of the runs, in which the 'control' data is normalized to 1 and all other data and their errors accordingly.
The following error types can be displayed:
\code{prop}	the propagated error (standard deviation);
\code{propSim}	the propagated error averaged from Monte Carlo simulation;
\code{eval}	the (unpropagated) error (standard deviation) resulting from averaging the expression evaluation;
\code{evalSim}	the (unpropagated) error (standard deviation) from the evaluated expressions of the Monte Carlo simulation.

If \code{type} = "bar", a normal error bar plot is displayed. If \code{type} = "hbar", ratios < 1 are displayed as hanging bars with (-1/value).
The latter is more appealing and visually discriminates the up-/down-regulated samples.
}

\value{
  A barplot with error bars displaying the ratios obtained from \code{\link{ratiocalc}}.  
}

\author{
  Andrej-Nikolai Spiess
}


\examples{
## quick (but not dirty!) analysis of qPCR data!
## normal bar plot with standard error of mean
## and t-test on crossing points
DAT <- modlist(reps, fluo = c(2:9), fct = l5())
GROUP <- c(1, 1, 2, 2, 3, 3, 4, 4)
res <-  ratiocalc(DAT, group = GROUP, ttest = "cp")
ratioplot(res, cex.names = 0.3, sem = TRUE, plot.t = "values")

## hanging bar plot with standard deviation
## and significance codes
ratioplot(res, cex.names = 0.3, type = "hbar", plot.t = "stars")

## same using the 'unpropagated' error and
## different y-scale
## (from averaged expression evaluations)
ratioplot(res, cex.names = 0.3, type = "hbar", 
	  errbar = "eval", ylim = c(-20, 59))

}

\keyword{nonlinear}

