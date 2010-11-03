\name{ratiobatch}
\alias{ratiobatch}

\title{Calculation of ratios in a batch format for multiple genes/samples}

\description{For multiple qPCR data from type 'pcrbatch', this function calculates ratios between samples, 
 using normalization against one or more reference gene(s), if supplied. The input can be single qPCR data or (more likely) data containing replicates. 
 This is essentially a version of \code{\link{ratiocalc}} that can handle multiple reference and control genes and multiple samples with replicates as found
 in large-scale qPCR runs such as 96- or 384-Well plates. The results are automatically stored as a file or copied into the clipboard.
}

\usage{
ratiobatch(data, group = NULL, plot = TRUE, 
           combs = c("same", "across", "all"), 
           type.eff = "mean.single", which.cp = "cpD2", 
           which.eff = "sli", dataout = "clip", ...)
}

\arguments{
  \item{data}{multiple qPCR data generated by \code{\link{pcrbatch}}.}
  \item{group}{a character vector defining the replicates (if any) as well as target and reference data. If \code{NULL}, it is taken from the column names of \code{data}. See 'Details'}.
  \item{plot}{logical. If \code{TRUE}, plots are displayed for the diagnostics and analysis.}
  \item{combs}{type of combinations between different samples (i.e. r1s1:g1s2). See 'Details'.}
  \item{type.eff}{type of efficiency averaging used. See \code{\link{ratiocalc}} for details.}
  \item{which.eff}{efficiency obtained from which method. See \code{\link{efficiency}} for details.}
  \item{which.cp}{threshold cycle obtained from which method. See \code{\link{efficiency}} for details.}   
  \item{dataout}{where to store the result dataframe. Either it is copied to the clipboard (default) or any path + filename such as \code{c:\\temp\\result.txt}.}   
  \item{...}{other parameters to be passed to \code{\link{ratiocalc}}.}
}

\details{
Similar to \code{\link{ratiocalc}}, the replicates of the 'pcrbatch' data columns are to be defined as a character vector with the following abbreviations:\cr

"g1s1":   gene-of-interest #1 in target sample #1\cr
"g1c1":   gene-of-interest #1 in control sample #1\cr
"r1s1":   reference gene #1 in target sample #1\cr
"r1c1":   reference gene #1 in control sample #1\cr

There is no distinction between the different technical replicates so that three different runs of gene of interest #1 in
 target sample #2 are defined as c("g1s2", "g1s2", "g1s2"). 

Examples:\cr
1 control with 2 gene-of-interests (2 technical replicates), 2 samples with 2 gene-of-interest (2 technical replicates):\cr
"g1c1", "g1c1", "g2c1", "g2c1", "g1s1", "g1s1", "g1s2", "g1s2", "g2s1", "g2s1", "g2s2", "g2s2"\cr
See also 'Examples'.  

The ratios are calculated for all pairwise 'rc:gc' and 'rs:gs' combinations according to:\cr
For all control samples \eqn{i = 1 \ldots I} and treatment samples \eqn{j = 1 \ldots J}, reference genes \eqn{k = 1 \ldots K} and gene-of-interests \eqn{l = 1 \ldots L}
calculate\cr

Without reference PCR: \deqn{\frac{E(g_is_j)^{cp(g_is_j)}}{E(g_kc_l)^{cp(g_kc_l)}}}
With reference PCR: \deqn{\frac{E(g_is_j)^{cp(g_is_j)}}{E(g_kc_l)^{cp(g_kc_l)}}/\frac{E(r_ms_n)^{cp(r_ms_n)}}{E(r_oc_p)^{cp(r_oc_p)}}}  

There are three different combination setups possible when calculating the pairwise ratios:

\code{combs = "same"}: reference genes, genes-of-interest, control and treatment samples are the \code{same}, i.e. 
 \eqn{i = k, m = o, j = n, l = p}.\cr
\code{combs = "across"}: control and treatment samples are the same, while the genes are combinated, i.e. \eqn{i \neq k, m \neq o, j = n, l = p, }.\cr
\code{combs = "all"}: reference genes, genes-of-interest, control and treatment samples are all combinated, i.e. 
 \eqn{i \neq k, m \neq o, j \neq n, l \neq p}.

The last setting rarely makes sense and is very time-intensive. \code{combs = "same"} is the most common setting, but \code{combs = "across"}
 also makes sense if different gene-of-interests and reference gene combinations should be calculated for the same samples.

Efficiencies can be taken from the individual curves or averaged from the replicates as described in the documentation to \code{\link{ratiocalc}}.
The different combinations of \code{type.eff}, \code{which.eff} and \code{which.cp} can yield very different results in ratio 
 calculation. We observed a relatively stable setup which minimizes the overall variance using the combination
  
\code{type.eff = "mean.single"}     # averaging efficiency across replicates\cr
\code{which.eff = "sli"}            # taking efficiency from the sliding window method\cr
\code{which.cp = "sig"}             # using the second derivative maximum of a sigmoidal fit\cr 

This is also the default setup in the function. The lowest variance can be obtained for the threshold cycles if the asymmetric 5-parameter \code{l5} model is used
 in the \code{\link{pcrbatch}} function.
}

\value{
A list with the following components:\cr
\item{resList}{a list with the results from the combinations as list items.}
\item{resDat}{a dataframe with the results in colums.}
Both \code{resList} and \code{resDat} have as names the combinations used for the ratio calculation.
If \code{plot = TRUE}, a plot matrix with the calculated mean ratios obtained from the Monte-Carlo simulation, permutation and error propagation is given as barcharts.
}

\author{
Andrej-Nikolai Spiess
}

\note{
This function can be used quite conveniently when the raw fluorescence data from the 96- or 384-well runs come from Excel with
 'Cycles' in the first column and run descriptions as explained above in the remaining column descriptions (such as 'r1c6'). Examples for
 a proper format can be found under \url{http://www.dr-spiess.de//qpcR//datasets.html}. This data may then be imported into \emph{R} by
 \code{dat <- pcrimport()}.
}

\examples{
\dontrun{
## one control, two samples, one gene-of-interest, 
## two reference genes, two replicates each
## replicates are averaged
DAT <- pcrbatch(reps, 2:19, l5)
GROUP <- c("r1c1", "r1c1", "r2c1", "r2c1", "g1c1", "g1c1",
           "r1s1", "r1s1", "r1s2", "r1s2", "r2s1", "r2s1",
           "r2s2", "r2s2", "g1s1", "g1s1", "g1s2", "g1s2") 
res <- ratiobatch(DAT, GROUP)    
## two controls, one sample, one gene-of-interest,
## one reference gene, no replicates
## use same efficiency E = 2         
DAT2 <- pcrbatch(reps, 2:7, l5)
GROUP2 <- c("r1c1", "r1c2", "g1c1", "g1c2", 
            "r1s1", "g1s1") 
res2 <- ratiobatch(DAT2, GROUP2, which.eff = 2) 
                   
## one control, one sample, three gene-of-interests,
## no reference gene, three replicates,
## using efficiency from sigmoidal model 
DAT3 <- pcrbatch(reps, 2:19, l5)
GROUP3 <- c("g1c1", "g1c1", "g1c1", "g2c1", "g2c1", "g2c1", "g3c1", "g3c1", "g3c1",
            "g1s1", "g1s1", "g1s1", "g2s1", "g2s1", "g2s1", "g3s1", "g3s1", "g3s1")
res3 <- ratiobatch(DAT3, GROUP3, which.eff = "sig")
}
}



\keyword{nonlinear}

