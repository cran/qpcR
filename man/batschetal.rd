\name{batschetal}
\alias{batsch1}
\alias{batsch2}
\alias{batsch3}
\alias{batsch4}
\alias{batsch5}

\title{qPCR dilution experiments from Batsch et al. (2008)}

\description{
High quality 4-fold dilution experiments with 5 dilution steps and 3 replicates each.
See 'Details' for the different setups.
}

\usage{
batsch1
batsch2
batsch3
batsch4
batsch5
}

\format{
Data frames with the PCR cycles and 15 qPCR runs with 3 replicates of five 4-fold dilutions.
The replicates are defined by FX.Y (X = dilution number, Y = replicate number).  
}

\details{
The real-time PCR was conducted with a Lightcycler 1.0 instrument using the following setups:\cr\cr
batsch1: Primers for rat SLC6A14, Taqman probes\cr
batsch2: Primers for human SLC22A13, Taqman probes\cr
batsch3: Primers for pig EMT, Taqman probes\cr
batsch4: Primers for chicken ETT, SybrGreen\cr 
batsch5: Primers for human GAPDH, SybrGreen\cr
}

\source{
Additional File 5 to the paper.
}

\references{
Simultaneous fitting of real-time PCR data with efficiency of amplification modeled as Gaussian function of target fluorescence.\cr
Batsch A et al., \emph{BMC Bioinformatics}, 2008, \bold{9}: 95.
}

\examples{
ml1 <- modlist(batsch1, model = l4)
plot(ml1)    
}

\keyword{models}
\keyword{nonlinear}
