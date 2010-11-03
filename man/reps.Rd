\name{reps}
\alias{reps}

\title{qPCR dilution experiments Set \#1 from the package author}

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
m1 <- pcrfit(reps, 1, 2, l5)
plot(m1)    
}

\keyword{models}
\keyword{nonlinear}
