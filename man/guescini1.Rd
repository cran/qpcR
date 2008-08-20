\name{guescini1}
\alias{guescini1}

\title{A qPCR dilution experiment with replicates (Lightcycler 480) from Guescini et al}

\description{
A high quality 10-fold dilution experiment with 7 dilution steps and 12 replicates each.
}

\usage{
data(guescini1)
}

\format{
A data frame with the PCR cycles and 84 qPCR runs with 12 replicates of seven 10-fold dilutions.
The replicates are defined by FX.Y (X = dilution number, Y = replicate number).  
}

\details{
The real-time PCR was conducted with primers for the NADH dehydrogenase 1 in a Lightcycler 480 (Roche).
The data is background subtracted.
}

\source{
Supplemental data 1 to the paper.
}

\references{
 A new real-time PCR method to overcome significant quantitative inaccuracy due to slight amplification inhibition.\cr
 Guescini M et al, \emph{BMC Bioinformatics}, 2008, \bold{9}: 326.
}

\examples{
\dontrun{
data(guescini1)
ml <- modlist(guescini1, fct = l5())
pcrplot(ml[[1]])
sapply(2:84, function(x) pcrplot(ml[[x]], add = TRUE))
}
}

\keyword{models}
\keyword{nonlinear}
