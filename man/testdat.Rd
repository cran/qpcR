\name{testdat}
\alias{testdat}

\title{dataset with noisy data to provoke fitting errors}

\description{
A dilution experiment with six 10-fold dilutions of the cDNA, and four replicates for each dilution. Each third and fourth replicate are made of noisy data that result in failed fits or correct fits but lack of sigmoidal structure, respectively. This dataset is mainly for testing the package functions for fail-safety against fitting errors and alike.
}

\usage{
testdat
}

\format{
A data frame with the PCR cycles and 24 qPCR runs with four replicates of seven 10-fold dilutions. The replicates are defined by FX.Y (X = dilution number, Y = replicate number). Each FX.3 has noisy data which fails to fit with the \code{\link{l5}} model, each FX.4 passes fitting but fails in sigmoidal structure detection by \code{\link{KOD}}.
}

\details{
The real-time PCR was conducted with primers for the S27a housekeeping gene in a Lightcycler 1.0 instrument (Roche Diagnostics).
}

\source{
Andrej-Nikolai Spiess & Nadine Mueller, Institute for Hormone and Fertlity Research, Hamburg, Germany.
}

\examples{
\dontrun{
## remove only non-fitted runs
ml1 <- modlist(testdat, model = l5, remove = "fit")
## remove sigmoidal outliers
ml2 <- modlist(testdat, model = l5, check = "uni2", remove = "KOD")  
}
}

\keyword{models}
\keyword{nonlinear}
