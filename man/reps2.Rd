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
data(reps2)### on single data
m1 <- multdrc(F1.1.1 ~ Cycles, data = reps2, fct = l3())
pcrplot(m1)
### on replicated data
repData <- repform(reps2[,1:4], c(0,1,1,1))
m2 <- multdrc(values ~ Cycles, curve = Curve, data = repData, fct = l5())
pcrplot(m2)
}

\keyword{models}
\keyword{nonlinear}
