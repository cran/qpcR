\name{REST}
\alias{REST}

\title{qpcR implementation of the popular REST software}

\description{
This function is in principle a port of the REST software to R (with minor modifications).
The paradigm behind the REST approach is to permutate threshold cycle values between controls and samples,
 calculate ratios for each permutation and finally construct a confidence interval on the obtained cumulative
 distribution. This approach has been replicated here with modifications as described under 'Details'.
}

\usage{
REST(eff.goi, eff.ref, data, B = 10000, alpha = 0.05, replace = FALSE, 
     bootstrap = NULL, perm = c("pairwise", "single"))
}

\arguments{
  \item{eff.goi}{the PCR efficiency of the gene-of-interest.}
  \item{eff.ref}{the PCR efficiency of the reference gene.}
  \item{data}{a dataframe with the threshold cycles for control ref, control goi, sample ref and sample goi. STRICTLY in that order.}  
  \item{B}{The number of permutations.}  	
  \item{alpha}{the confidence level.}
  \item{replace}{logical. Should permutations be conducted with replacement?}
  \item{bootstrap}{The proportion of data to be used for the permutation (bootstrapping). Must be between 0 and 1, default is 1.}
  \item{perm}{type of permutation. 'pairwise' (default) is the usual REST procedure ('pairwise fixed reallocation'). See 'Details'.}
}

\details{
The following modifications were made to (enhance) the original REST version.\cr
a)  Permutation sampling can be controlled with/without replacement and data subsetting (bootstrap).\cr
b)  Additionally to the pairwise fixed reallocation, a permutation regime without fixing control and GOI of reference/sample
    can be conducted. Thus, all four input columns of \code{data} are permutated. This increases the number of non-redundant
    permutations, but contradicts the pradigm of a pairwise observation (two PCRs on the same sample). The user may judge, which
    approach is more sensible.\cr
c)  The list item \code{uniques} reveals the redundancy of permutations (=> no information gain for confidence interval).\cr
}

\value{
A list with the following components:
\item{ratios}{the sorted ratios obtained from the permutations.}
\item{conf}{the confidence interval for \code{alpha}.}
\item{uniques}{a dataframe with the UNIQUE permutations. Reveals redundancy in dependence of \code{B}.}           
}

\author{
Andrej-Nikolai Spiess
}

\references{
Relative Expression Software Tool (REST) for group wise comparison and
statistical analysis of relative expression results in real-time PCR.\cr
Pfaffl MW, Horgan GW & Dempfle L.\cr
\emph{Nucl Acids Res} 2002, \bold{30}: E36.
}

\examples{
### artificial data with very low variance 
### in threshold cycles (~0.5%)
### See how this still gives a 
### relatively large confidence interval!
DATA <- data.frame(ref.C = rnorm(10, 26 , 0.1),
                   goi.C = rnorm(10, 27, 0.1),
                   ref.S = rnorm(10, 26, 0.1),
                   goi.S = rnorm(10, 24, 0.1))

res <- REST(eff.goi = 2.00, eff.ref = 1.97, data = DATA, 
            B = 10000, perm = "pairwise") 
            
boxplot(res$ratios, col = 2, outline = FALSE)
}

\keyword{models}
\keyword{nonlinear}
