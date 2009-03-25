\name{pcrfit}
\alias{pcrfit}

\encoding{latin1}

\title{Simple wrapper function for fast model fitting}

\description{
This is a wrapper function with easy syntax which simplifies the model fitting process by \code{\link[drc]{multdrc}}.
}

\usage{
pcrfit(data, cyc = 1, fluo, fct = l5())
}

\arguments{
  \item{data}{the name of the dataframe containing the variables.}
  \item{cyc}{the column number with the cycle data. Defaults to 1.}
  \item{fluo}{the column containing the raw fluorescence data.}
  \item{fct}{the model to be used for the analysis, with the usual 'drc' package nomenclature. Defaults to 'l5'.}  
}

\details{
This function is to be used at the single run level. Otherwise use \code{\link{pcrbatch}}.
}

\value{
  A 'drc' model.
} 

\author{
  Andrej-Nikolai Spiess
}


\references{
Ritz C & Streibig JC. Bioassay analysis using R. \emph{J Stat Soft} (2005), \bold{12}: 1-22.
}

\seealso{
  The function \code{\link[drc]{multdrc}} for more detailed model specification.
}

\examples{
### simple l4() fit of F1.1 of the 'reps' dataset
m <- pcrfit(reps, 1, 2, l4())
pcrplot(m)

### same with five-parameter model
m2 <- pcrfit(reps, 1, 2)
pcrplot(m2, add = TRUE, col = 2)
}

\keyword{models}
\keyword{nonlinear}
