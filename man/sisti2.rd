\name{sisti2}
\alias{sisti2}

\title{qPCR dilution experiments Set #2 from Sisti et al. (2010)}

\description{
A high quality 2-fold inhibitor dilution experiment with 9 dilution steps of tannic acid and 6 replicates each.
}

\usage{
data(sisti2)
}

\format{
A data frame with the PCR cycles and 54 qPCR runs with 6 replicates of nine 2-fold dilutions of tannic acid.
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
data(sisti2)
ml <- modlist(sisti2, model = l4)
plot(ml) 
}
}

\keyword{models}
\keyword{nonlinear}
