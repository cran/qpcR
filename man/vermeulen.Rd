\name{vermeulenetal}
\alias{vermeulen1}
\alias{vermeulen2}
\encoding{latin1}

\title{A large-scale qPCR dataset from Vermeulen et al. (2009)}

\description{
A subset of a larger dataset from Vermeulen \emph{et al}. (2009). The original dataset consists of 64 genes (59 prognostic genes and 5 reference genes) with 384 runs/gene. The 'qpcR' package includes a subset of the first 20 runs for each of the 64 genes, hence 1280 runs (\code{vermeulen1}) and the corresponding dilution data for all 64 genes (\code{vermeulen2}). 
}

\usage{
vermeulen1
vermeulen2
}

\format{
\code{vermeulen1}: A data frame with 1280 runs named by X.Y, with X = gene name and Y = sample number.\cr
\code{vermeulen2}: A data frame with dilution data for all 64 genes, with five 10-fold dilutions and 3 replicates each. Data is in the format X.STD_Y.Z, with X = gene name, Y = copy number and Z = replicate number. 
}

\details{
The real-time PCR was conducted in a Lightcycler 480 (Roche) using SybrGreen I chemistry. The investigated genes are AHCY, AKR1C1, ALUsq(Eurogentec), ARHGEF7, BIRC5, CAMTA1, CAMTA2, CD44, CDCA5, CDH5, CDKN3, CLSTN1, CPSG3, DDC, DPYSL3, ECEL1, ELAVL4, EPB41L3, EPHA5, EPN2, FYN, GNB1, HIVEP2, HMBS, HPRT1, IGSF4, INPP1, MAP2K4, MAP7, MAPT, MCM2, MRPL3, MTSS1, MYCN(4), NHLH2, NM23A, NRCAM, NTRK1, ODC1, PAICS, PDE4DIP, PIK3R1, PLAGL1, PLAT, PMP22, PRAME, PRDM2, PRKACB, PRKCZ, PTN, PTPRF, PTPRH, PTPRN2, QPCT, SCG2, SDHA(1), SLC25A5, SLC6A8, SNAPC1, TNFRSF, TYMS, UBC(2), ULK2 and WSB1.
}

\source{
Originally, raw data was available at \url{http://medgen.ugent.be/jvermeulen}, but site is down. The complete (\code{vermeulen_all}) and smaller (\code{vermeulen_sub}) datasets can be downloaded from \url{http://www.dr-spiess.de/qpcR/datasets.html}. 
}

\references{
Predicting outcomes for children with neuroblastoma using a multigene-expression signature: a retrospective SIOPEN/COG/GPOH study.\cr
Vermeulen J, De Preter K, Naranjo A, Vercruysse L, Van Roy N, Hellemans J,
Swerts K, Bravo S, Scaruffi P, Tonini GP, De Bernardi B, Noguera R, Piqueras M,
Cañete A, Castel V, Janoueix-Lerosey I, Delattre O, Schleiermacher G, Michon J,
Combaret V, Fischer M, Oberthuer A, Ambros PF, Beiske K, Bénard J, Marques B,
Rubie H, Kohler J, Pötschger U, Ladenstein R, Hogarty MD, McGrady P, London WB,
Laureys G, Speleman F & Vandesompele J.\cr 
\emph{Lancet Oncol} (2009), \bold{10}:663-671.
}

\examples{
\dontrun{
## make model list and plot
## with colors for each gene 
COL <- rep(rainbow(64), each = 20)
ml <- modlist(vermeulen, model = l4)
plot(ml, col = COL)
}
}

\keyword{models}
\keyword{nonlinear}
