\name{S27}
\alias{S27}

\title{PCR data from 15 different testicular biopsies (Roche Lightcycler)}

\description{
The data was obtained from cDNA reverse transcribed from the isolated RNA of 15 different testicular biopsies.
}

\usage{
data(S27)
}

\format{
A data frame with cycle number and 15 different qPCR runs (F1 - F15).  
}

\details{
The real-time PCR was conducted with primers for the S27a housekeeping gene in a Lightcycler 1.0 
instrument (Roche Diagnostics).
}

\source{
Andrej-Nikolai Spiess & Caroline Feig, Department of Andrology, University Hospital Hamburg, Germany.
}

\examples{
data(S27)
### on single data
m1 <- multdrc(F1 ~ Cycles, data = S27, fct = l5())
pcrplot(m1)
### defining 15 replicates
repData <- repform(S27[,1:16], c(0, rep(1, 15)))
m2 <- multdrc(values ~ Cycles, curve = Curve, data = repData, fct = l5())
pcrplot(m2, type = "all")
}

\keyword{models}
\keyword{nonlinear}
