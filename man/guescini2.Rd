\name{guescini2}
\alias{guescini2}

\title{qPCR dilution experiments Set \#2 from Guescini et al. (2008)}

\description{
A high quality experiment in which a decreasing amount of PCR mix mimics PCR inhibition.
}

\usage{
data(guescini2)
}

\format{
A data frame with the PCR cycles and 60 qPCR runs from 12 replicates with 5 decreasing steps of PCR mix.
The replicates are defined by FX.Y (X = PCR mix dilution number, Y = replicate number). 
For example:\cr\cr
F1.X	100\% Mix\cr
F2.X	90\%  Mix\cr
F3.X	80\%  Mix\cr
F4.X 	70\%  Mix\cr
F5.X	60\%  Mix\cr
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
## effect of decreasing mix
## on PCR efficiency
data(guescini2)
ml <- modlist(guescini2, model = l4)
effs <- sapply(ml, function(x) efficiency(x, plot = FALSE)$eff)

## mean for the replicates
tapply(effs, gl(5, 12), function(x) mean(x, na.rm = TRUE))
}

\keyword{models}
\keyword{nonlinear}
