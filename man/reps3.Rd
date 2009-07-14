\name{reps3}
\alias{reps3}

\title{A qPCR dilution experiment with replicates (Stratagene MXPro3000P)}

\description{
A high quality dilution experiment with six 4-fold dilutions of the cDNA, and three replicates for each dilution.
}

\usage{
data(reps3)
}

\format{
A data frame with the PCR cycles and 21 qPCR runs with three replicates of six 4-fold dilutions.
The replicates are defined by FX.Y (X = dilution number, Y = replicate number).  
}

\details{
The real-time PCR was conducted with primers for the S27a housekeeping gene in a MXPro3000P instrument (Stratagene).
Data was ROX-normalized, but without smoothing.
}

\source{
Heike Cappallo-Obermann, Bone Marrow Transplantation Unit, University Hospital Hamburg-Eppendorf.
}

\examples{
## on single data
data(reps3)
m1 <- pcrfit(reps3, 1, 2, l5)
pcrplot(m1)  
}

\keyword{models}
\keyword{nonlinear}
