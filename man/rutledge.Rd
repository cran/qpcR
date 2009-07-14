\name{rutledge}
\alias{rutledge}

\title{qPCR dilution experiment with replicates and different runs from Rutledge et al (Biorad Opticon)}

\description{
A dilution experiment of a 102 bp amplicon with six 10-fold dilutions, with four replicates each in five independent runs.
}

\usage{
data(rutledge)
}

\format{
A data frame with the PCR cycles and and the six 10-fold dilutions, with four replicates each and five independent runs,
organized as X.Ry.Z, with x = the dilution number, y = the run number and z = the replicate number.
}

\details{
The real-time PCR was conducted for a 102 bp amplicon in an Opticon2 instrument.
}

\source{
Supplemental data 1 to the paper.
}

\references{
 Sigmoidal curve-fitting redefines quantitative real-time PCR with the prospective of developing automated high-throughput applications.
 Rutledge RG, \emph{Nucleic Acids Research}, 2004, \bold{32}, e178.
}

\examples{
data(rutledge)
m1 <- pcrfit(rutledge, 1 , 2, l5)
pcrplot(m1)
}

\keyword{models}
\keyword{nonlinear}
