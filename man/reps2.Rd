\name{reps2}
\alias{reps2}

\title{Another qPCR dilution experiment with replicates (Roche Lightcycler)}

\description{
A dilution experiment of two different cDNAs with five 4-fold dilutions of the cDNA, and three replicates for each dilution.
}

\usage{
data(reps2)
}

\format{
A data frame with the PCR cycles and 30 qPCR runs with three replicates of five 4-fold dilutions.
The replicates are defined by FX.Y.Z (X = cDNA number, Y = dilution number, Z = replicate number).  
}

\details{
The real-time PCR was conducted with primers for the S27a housekeeping gene in a Lightcycler 1.0 
instrument (Roche Diagnostics).
}

\source{
Heike Cappallo-Obermann, Bone Marrow Transplantation Unit, University Hospital Hamburg-Eppendorf.
}

\examples{ 
data(reps2)
m1 <- pcrfit(reps2, 1, 2, l3)
plot(m1) 
}

\keyword{models}
\keyword{nonlinear}
