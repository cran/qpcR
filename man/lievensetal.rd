\name{lievensetal}
\alias{lievens1}
\alias{lievens2}
\alias{lievens3}

\title{qPCR dilution and inhibition data from Lievens et al. (2012)}

\description{
\code{lievens1}: High quality 5-fold dilution experiments with 5 dilution steps and 18 replicates.\cr 
\code{lievens2}: Inhibition data with five different concentrations of isopropanol and 18 replicates.\cr
\code{lievens3}: Inhibition data with five different amounts of tannic acid per reaction and 18 replicates.\cr
}

\usage{
lievens1
lievens2
lievens3
}

\format{
\code{lievens1}: 90 qPCR runs with five 5-fold dilutions and 18 replicates each. Named by SX.Y, with X = dilution step and Y = replicate.\cr
\code{lievens2}: 90 qPCR runs with five different concentrations of isopropanol (2.5\%, 0.5\%, 0.1\%, 0.02\% and 0.004\% (v/v)) and 18 replicates each. Named by SX.Y, with X = concentration step and Y = replicate.\cr
\code{lievens3}: 90 qPCR runs with five different amounts of tannic acid per reaction (5 ng, 1 ng, 0.2 ng, 0.04 ng and 0.008 ng) and 18 replicates each. Named by SX.Y, with X = concentration step and Y = replicate.\cr
}

\details{
The real-time PCR was conducted with an ABI7300 (ABI) or Biorad IQ5 (Biorad) instrument using SybrGreen I chemistry and primers for the soybean lectin endogene Le1.
}

\source{
Supplementary Data to the paper.
}

\references{
Enhanced analysis of real-time PCR data by using a variable efficiency model: FPK-PCR.\cr
Lievens A, Van Aelst S, Van den Bulcke M & Goetghebeur E.\cr
\emph{Nucleic Acids Res} (2012), \bold{40}:e10.
}

\examples{
\dontrun{
## lievens1
ml1 <- modlist(lievens1, model = l4)
plot(ml1) 

## lievens2
COL <- rep(1:4, each = 18)
ml2 <- modlist(lievens2, model = l4)
plot(ml2, col = COL)
}
}

\keyword{models}
\keyword{nonlinear}
