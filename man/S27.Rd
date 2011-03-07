\name{S27}
\alias{S27}

\title{qPCR data of the S27a transcript in 15 testicular biopsies}

\description{
The data was obtained from cDNA reverse transcribed from the isolated RNA of 15 different testicular biopsies.
}

\usage{
S27
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
m1 <- pcrfit(S27, 1, 2, l5)
plot(m1)     
}

\keyword{models}
\keyword{nonlinear}
