\name{LOF.test}
\alias{LOF.test}

\title{A formal Lack-Of-Fit test of a nonlinear model against a one-way ANOVA model}

\description{
Tests the nonlinear model against a more general one-way ANOVA model and from a likelihood ratio test.
P-values are derived from the F- and \eqn{\chi^2} distribution, respectively. 
}

\usage{
LOF.test(object)
}

\arguments{
  \item{object}{an object of class 'replist' or 'nls', which was fit with replicate response values.}
}

\value{
A list with the following components:
  \item{pF}{the p-value from the F-test against the one-way ANOVA model.}
  \item{pLR}{the p-value from the likelihood ratio test against the one-way ANOVA model.}      
}

\details{
The one-way ANOVA model is constructed from the \code{data} component of the nonlinear model by factorizing each of the predictor values.
Hence, the nonlinear model becomes a submodel of the one-way ANOVA model and we test both models with the null hypothesis that the ANOVA model
 can be simplified to the nonlinear model (Lack-of-fit test). This is done by two approaches:

1) an F-test (Bates & Watts, 1988).\cr
2) a likelihood ratio test (Huet et al., 2004).

P-values are derived from an F-distribution (1) and a \eqn{\chi^2} distribution (2).
}
                

\author{
Andrej-Nikolai Spiess
}    

\references{
Bates DM & Watts DG.\cr
Nonlinear Regression Analysis and its Applications.\cr
John Wiley & Sons (1988), New York.\cr

Huet S et al.\cr
Statistical Tools for Nonlinear Regression: A Practical Guide with S-PLUS and R Examples.\cr
Springer Verlag (2004), New York, 2nd Ed.
}  

\examples{
## Example with a 'replist'
## no lack-of-fit
ml <- modlist(reps, fluo = 2:5, model = l5)
rl <- replist(ml, group = c(1, 1, 1, 1))
LOF.test(rl)

## Example with a 'nls' fi
## there is a lack-of-fit
DNase1 <- subset(DNase, Run == 1)
fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1) 
LOF.test(fm1DNase1)
}

\keyword{models}
\keyword{nonlinear}
