\name{rutledge}
\alias{rutledge}

\title{qPCR dilution experiments from Rutledge (2004)}

\description{
A dilution experiment of a 102 bp amplicon with six 10-fold dilutions, with four replicates each in five independent runs.
}

\usage{
rutledge
}

\format{
A data frame with the PCR cycles and and the six 10-fold dilutions, with four replicates each and five independent runs, organized as X.Ry.Z, with X = the dilution number, Y = the run number and Z = the replicate number.
}

\details{
The real-time PCR was conducted for a 102 bp amplicon in an Opticon2 instrument.
}

\source{
Supplemental data 1 to the paper.
}

\references{
Sigmoidal curve-fitting redefines quantitative real-time PCR with the prospective of developing automated high-throughput applications.\cr
Rutledge RG.\cr
\emph{Nucleic Acids Research} (2004), \bold{32}: e178.
}

\examples{
data(rutledge)
m1 <- pcrfit(rutledge, 1 , 2, l5)
plot(m1)
}

\keyword{models}
\keyword{nonlinear}
