\name{guescinietal}
\alias{guescini1}
\alias{guescini2}


\title{qPCR dilution experiments from Guescini et al. (2008)}

\description{
\code{guescini1}: A high quality 10-fold dilution experiment with 7 dilution steps and 12 replicates each.\cr
\code{guescini2}: A high quality experiment in which a decreasing amount of PCR mix mimics PCR inhibition.
}

\usage{
guescini1
guescini2
}

\format{
\code{guescini1}: A data frame with the PCR cycles and 84 qPCR runs with 12 replicates of seven 10-fold dilutions. The replicates are defined by FX.Y (X = dilution number, Y = replicate number).\cr
\code{guescini2}: A data frame with the PCR cycles and 60 qPCR runs from 12 replicates with 5 decreasing steps of PCR mix. The replicates are defined by FX.Y (X = PCR mix dilution number, Y = replicate number).
}

\details{
The real-time PCR was conducted with primers for the NADH dehydrogenase 1 in a Lightcycler 480 (Roche). The data is background subtracted.
}

\source{
Supplemental data 1 to the paper.
}

\references{
A new real-time PCR method to overcome significant quantitative inaccuracy due to slight amplification inhibition.\cr
Guescini M, Sisti D, Rocchi MB, Stocchi L & Stocchi V.\cr
\emph{BMC Bioinformatics} (2008), \bold{9}: 326.\cr
}

\examples{
\dontrun{
## dilution set with replicates
ml1 <- modlist(guescini1, model = l4)
plot(ml1) 

## effect of decreasing mix
## on PCR efficiency,
## mean for the replicates
ml2 <- modlist(guescini2, model = l4)
effs <- sapply(ml2, function(x) efficiency(x, plot = FALSE)$eff)
tapply(effs, gl(5, 12), function(x) mean(x, na.rm = TRUE))
}
}

\keyword{models}
\keyword{nonlinear}
