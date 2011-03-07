\name{getPar}
\alias{getPar}

\title{Batch calculation of qPCR efficiencies/threshold cycles with simple output, especially tailored to high-throughput data}

\description{
This is a cut-down version of \code{\link{pcrbatch}}, starting with data of class 'modlist', which delivers a simple two-row dataframe output, with the calculated threshold cycles in the first row
 and the efficiencies in the second row. The column names are deduced from the run names. All calculations have been error-protected through \code{\link{try}}, so whenever there is any kind of error
 (fitting, efficiency estimation etc), \code{NA} is returned. This function can be used with high throughput data quite conveniently. All methods as in
 \code{\link{pcrbatch}} are available. The results are automatically copied to the clipboard.  
}

\usage{
getPar(x, cp = "cpD2", eff = "sigfit", ...)
}

\arguments{
  \item{x}{an object of class 'pcrfit' or 'modlist'.}
  \item{cp}{which method for threshold cycle estimation. Any of the methods in \code{\link{efficiency}}, i.e. "cpD2" (default), "cpD1", "maxE", "expR", "Cy0", "CQ", "maxRatio".}
  \item{eff}{which method for efficiency estimation. Either "sigfit" (default), "sliwin" or "expfit".}  
  \item{...}{other parameters to be passed to \code{\link{efficiency}}, \code{\link{sliwin}} or \code{\link{expfit}}.} 
}

\details{
Takes about 4 sec for 100 runs on a Pentium 4 Quad-Core (3 Ghz), INCLUDING the 'modlist' creation.
}

\value{
A dataframe, which is automatically copied to the clipboard, with threshold cycle values in the first row and efficiencies in the second.
}

\author{
Andrej-Nikolai Spiess.
}

\examples{
## simple example with
## plotting of threshold cycles
\dontrun{
ml <- modlist(rutledge, model = l5)
res <- getPar(ml, cp = "cpD2", eff = "sliwin")
barplot(res[1, ], las = 2)
}
}

\keyword{models}
\keyword{nonlinear}
