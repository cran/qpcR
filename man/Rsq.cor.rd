\name{Rsq.cor}
\alias{Rsq.cor}

\title{R-square derived from correlation of data and fitted values}

\description{
Calculates the R-square value for objects of class \code{nls}, \code{lm}, \code{glm}, \code{drc}
 or any other models from which \code{\link{fitted}} and \code{\link{residuals}} can be extracted.
The value is calculated as a correlation measure. See 'Details'.
}

\usage{
Rsq.cor(object)
}

\arguments{
  \item{object}{a fitted model.}
 }

\value{
The R-square value of the fit.
}

\details{
Calculates the R-square by \deqn{R_{cor}^2 = cor(y_i, \hat{y}_i)^2}
}

\author{
Andrej-Nikolai Spiess
}

\references{
Summarizing the predictive power of a generalized linear model.\cr
B. Zheng & A. Agresti\cr
\emph{Statist. Med.} (2000), \bold{19}: 1771 - 1781.
}



\examples{
m <- pcrfit(reps, 1, 2, l5)
Rsq.cor(m)

## compare to 'standard' R-square
Rsq(m) 
}

\keyword{models}
\keyword{nonlinear}
