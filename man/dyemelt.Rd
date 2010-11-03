\name{dyemelt}
\alias{dyemelt}

\encoding{latin1}

\title{Melting curves of a 4-plex qPCR with different fluorescence dyes}

\description{
A melting curve analysis of a multiplex real-time PCR on genomic DNA, producing 4 amplicons.
The fluorescence dyes used were EvaGreen, SybrGreen I and Syto-13.
}

\usage{
data(dyemelt)
}

\format{
A data frame with the temperature values in columns 1, 2, 3 and the corresponding fluorescence data in columns 2, 4, 6. 
}

\details{
The melting curve was conducted with AZF deletion-specific primers in a Lightcycler 1.0 instrument (Roche).
}

\examples{
data(dyemelt)
plot(dyemelt[, 1], dyemelt[, 2], xlab = "Temperature [°C]", 
     ylab = "Raw fluorescence")  
}

\keyword{models}
\keyword{nonlinear}     