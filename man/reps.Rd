\name{reps}
\alias{reps}

\title{qPCR dilution experiment with replicates (Roche Lightcycler)}

\description{
A dilution experiment with seven 10-fold dilutions of the cDNA, and four replicates for each dilution.
}

\usage{
data(reps)
}

\format{
A data frame with the PCR cycles and 28 qPCR runs with four replicates of seven 10-fold dilutions.
The replicates are defined by F1.1 - F1.4 (first dilution), F2.1 - F2.4 (second diltution) etc.  
}

\details{
The real-time PCR was conducted with primers for the S27a housekeeping gene in a Lightcycler 1.0 
instrument (Roche Diagnostics).
}

\source{
Andrej-Nikolai Spiess & Nadine Mueller, Institute for Hormone and Fertlity Research, Hamburg, Germany.
}

\examples{
data(reps)
### on single data
m1 <- pcrfit(reps, 1, 2, l3())
pcrplot(m1)
### on replicated data
repData <- repform(reps[,1:5], c(0,1,1,1,1))
m2 <- drmfit(values ~ Cycles, curve = Curve, data = repData, fct = l5())
pcrplot(m2)
}

\keyword{models}
\keyword{nonlinear}
