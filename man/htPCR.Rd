\name{htPCR}
\alias{htPCR}

\title{High throughput qPCR experiment with 8858 runs}

\description{
A high throughput experiment containing 8858 runs from a 95 x 96 PCR grid.
Some wells have been removed, probably due to failed amplification.
}

\usage{
data(htPCR)
}

\format{
A data frame with 8858 PCR runs.  
}

\details{
The real-time PCR was conducted with proprietary primers in a Fluidigm instrument (Biomark) using EvaGreen and ROX normalization.   
}

\source{
Kindly supplied by Roman Bruno.
}

\examples{
\dontrun{
data(htPCR)
ml1 <- modlist(htPCR, 1, 2:200, model = l5)
getPar(ml1, type = "curve") 
}
}

\keyword{models}
\keyword{nonlinear}
