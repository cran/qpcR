\name{boggy}
\alias{boggy}

\title{qPCR dilution experiments from Boggy et al. (2010)}

\description{
A dilution experiment with six 10-fold dilutions of a synthetic template, and two replicates for each dilution.
}

\usage{
data(boggy)
}

\format{
A data frame with the PCR cycles and 12 qPCR runs with two replicates of six 10-fold dilutions.
The replicates are defined by F1.1 - F1.2 (first dilution), F2.1 - F2.2 (second diltution) etc.  
}

\details{
The real-time PCR was conducted with primers for a synthetic template, consisting of a secondary structure-optimized
 random sequence (129 bp), in a Chromo4 instrument (BioRad) with Syto-13 dye.
}

\source{
Additional File S1 to the paper.
}

\references{
A Mechanistic Model of PCR for Accurate Quantification of Quantitative PCR Data.\cr
Boggy GJ and Woolf PJ.\cr
\emph{PLOS One}, 2010, \bold{5}: e12355.
}

\examples{
data(boggy)
m1 <- pcrfit(boggy, 1, 2, l5)
plot(m1)    
}

\keyword{models}
\keyword{nonlinear}
