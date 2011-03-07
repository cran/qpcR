\name{repsdat}
\alias{reps}
\alias{reps2}
\alias{reps3}

\title{qPCR dilution experiments from the package author}

\description{
\code{reps}: A dilution experiment with seven 10-fold dilutions of the cDNA, and four replicates for each dilution.\cr
\code{reps2}: A dilution experiment of two different cDNAs with five 4-fold dilutions of the cDNA, and three replicates for each dilution.\cr
\code{reps3}: A dilution experiment with six 4-fold dilutions of the cDNA, and three replicates for each dilution.
}

\usage{
reps
reps2
reps3
}

\format{
\code{reps}: A data frame with the PCR cycles and 28 qPCR runs with four replicates of seven 10-fold dilutions. The replicates are defined by FX.Y (X = dilution number, Y = replicate number).\cr
\code{reps2}: A data frame with the PCR cycles and 30 qPCR runs with three replicates of five 4-fold dilutions. The replicates are defined by FX.Y.Z (X = cDNA number, Y = dilution number, Z = replicate number).\cr
\code{reps3}: A data frame with the PCR cycles and 21 qPCR runs with three replicates of six 4-fold dilutions. The replicates are defined by FX.Y (X = dilution number, Y = replicate number).
}

\details{
The real-time PCR was conducted with primers for the S27a housekeeping gene in a Lightcycler 1.0 instrument (Roche Diagnostics) for \code{reps} and \code{reps2} or in a MXPro3000P instrument (Stratagene) for \code{reps3}. \code{reps3} was ROX-normalized.
}

\source{
Andrej-Nikolai Spiess & Nadine Mueller, Institute for Hormone and Fertlity Research, Hamburg, Germany.
}

\examples{
m1 <- pcrfit(reps, 1, 2, l5)
plot(m1)
}

\keyword{models}
\keyword{nonlinear}
