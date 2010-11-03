\name{sisti1}
\alias{sisti1}

\title{qPCR dilution experiments Set \#1 from Sisti et al. (2010)}

\description{
A high quality 10-fold dilution experiment with 6 dilution steps and 12 replicates each.
}

\usage{
data(sisti1)
}

\format{
A data frame with the PCR cycles and 72 qPCR runs with 12 replicates of six 10-fold dilutions.
The replicates are defined by FX.Y (X = dilution number, Y = replicate number).  
}

\details{
The real-time PCR was conducted with primers for the MT-ND1 in a Lightcycler 480 (Roche).
The data is background subtracted.
}

\source{
Supplemental data 1 to the paper.
}

\references{
Shape based kinetic outlier detection in real-time PCR.\cr
Sisti D et al, \emph{BMC Bioinformatics}, 2010, \bold{11}: 186. 
}

\examples{
\dontrun{
data(sisti1)
ml <- modlist(sisti1, model = l4)
plot(ml) 
}
}

\keyword{models}
\keyword{nonlinear}
