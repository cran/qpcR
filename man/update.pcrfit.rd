\name{update.pcrfit}
\alias{update.pcrfit}

\title{Updating and refitting a qPCR model}

\description{
Updates and re-fits a model of class 'pcrfit'. 
}

\usage{
\method{update}{pcrfit}(object, ..., evaluate = TRUE)
}

\arguments{
  \item{object}{a fitted model of class 'pcrfit'.}
  \item{...}{arguments to alter in object.}
  \item{evaluate}{logical. If \code{TRUE}, model is re-fit; otherwise an unevaluated call is returned.}   
}

\value{
An updated model of class 'pcrfit' and 'nls'.
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
The function \code{\link{pcrfit}} in this package.
}

\examples{
m <- pcrfit(reps, 1, 2, l4)

## update model
update(m, model = baro5)

## update qPCR run
update(m, fluo = 20)

## update data
update(m, data = guescini1) 

## update 'optim' method
update(m, opt.method = "GA")   
}

\keyword{models}
\keyword{nonlinear}
