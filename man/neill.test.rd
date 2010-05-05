\name{neill.test}
\alias{neill.test}
\encoding{latin1}

\title{Neill's lack-of-fit test when replicates are lacking}  

\description{
This is a test for nonlinear models in the absence of replicates. A \code{grouping} of the predictor values has to be provided. If missing, it is calculated from
 cutting the dendrogram into groups of at least 2 predictor values. See 'References' for details. 
 Works on any model WITHOUT replicates of class 'pcrfit' or 'nls'.
}

\usage{
neill.test(object, grouping) 
}

\arguments{
  \item{object}{an object of class 'pcrfit'.}
  \item{grouping}{vector that provides the grouping of the predictor ('Cycles') values.}  
}

\details{
Let the total number of observations be denoted by \eqn{\sum{n_i}} and let \eqn{\bar{G}(\theta)} be defined by 
\deqn{\bar{G}(\theta) = (G(x_1, \theta)j_{n_{1}}', \ldots,G(x_M, \theta)j_{n_{M}}')'},
where \eqn{j_{n_{i}}} is an \eqn{n_i \times 1} vector of ones, \eqn{i = 1,2, \ldots, M}. Also, \eqn{Y = (Y_{11}, \ldots, Y_{1n_{1}}, \ldots, Y_{Mn_{m}})'}
 and \eqn{\bar{Y} = (\bar{Y}_{1 \cdot j_{n_{1}}'}, \ldots, \bar{Y}_{M \cdot j_{n_{M}}'})'}.
 To test \deqn{H_0: E(Y) = \bar{G}(\theta)} vs. \deqn{H_{\alpha}: E(Y) \neq \bar{G}(\theta)} let the statistic \eqn{F} be defined by
 \deqn{F = [(N - M)/(M - p)](\parallel \bar{Y} - \bar{G}(\hat{\theta}) \parallel^2 / \parallel Y - \bar{Y} \parallel^2)},
 where \eqn{\hat{\theta}} is the least-squares parameter estimator of \eqn{\theta}. Reject \eqn{H_0} if observed \eqn{F > F_{M-p, N-M}^\alpha}.
 This is a nonlinear analogue to the Lack-of-fit test in linear models with replication.
}

\value{
The p-value from the test.
}

\author{
Andrej-Nikolai Spiess, taken in part from function \code{neill.test} of the 'drc' package.
}

\references{
## Original publication:
\cr
Neill JW.\cr
Testing for lack-of-fit in nonlinear regression.\cr
\emph{Ann Statist} (1988), \bold{16}, 733-740.\cr
## Using dendrogram splitting for grouping:
\cr
Ritz C, Martinussen T.\cr
Lack-of-fit tests for assessing mean structures for continuous dose-response data.\cr
\emph{Environ Ecol Stat} (2010), \bold{xx}, xxx-xxx.  
}

\examples{
## compare two models
m1 <- pcrfit(reps, 1, 2, l4)
m2 <- pcrfit(reps, 1, 2, l5)
neill.test(m1)
neill.test(m2)

\dontrun{
## using example from 'nls'
## Fails when replicates are given
DNase1 <- subset(DNase, Run == 1)  
fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
neill.test(fm1DNase1)
## But works if replicates are removed
DNase2 <- DNase1[-c(2, 4, 6, 8, 10, 12, 14, 16),]
fm1DNase2 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase2)
neill.test(fm1DNase2) 
} 
}

\keyword{models}
\keyword{nonlinear}
