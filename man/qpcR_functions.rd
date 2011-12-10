\name{qpcR_functions}
\alias{l7}
\alias{l6}
\alias{l5}
\alias{l4}
\alias{l3}
\alias{b7}
\alias{b6}
\alias{b5}
\alias{b4}
\alias{b3}   
\alias{expGrowth}  
\alias{mak2}
\alias{mak3}
\alias{mak3n}
\alias{chag}


\title{The nonlinear/mechanistic models implemented in qpcR}

\description{
A summary of all available models implemented in this package.
}

\usage{
l7
l6
l5
l4
l3
b7
b6
b5
b4
b3
expGrowth 
mak2
mak3
mak3n
chag
}

\details{
The following nonlinear sigmoidal models are implemented:\cr\cr
\bold{l7:} \deqn{f(x) = c + k1 \cdot x + k2 \cdot x^2 + \frac{d-c}{(1+exp(b(log(x)-log(e))))^f}}
\bold{l6:} \deqn{f(x) = c + k \cdot x + \frac{d-c}{(1+exp(b(log(x)-log(e))))^f}}
\bold{l5:} \deqn{f(x) = c + \frac{d-c}{(1+exp(b(log(x)-log(e))))^f}}
\bold{l4:} \deqn{f(x) = c + \frac{d-c}{1+exp(b(log(x)-log(e)))}}
\bold{l3:} \deqn{f(x) = \frac{d}{1+exp(b(log(x)-log(e)))}} 
\bold{b7:} \deqn{f(x) = c + k1 \cdot x + k2 \cdot x^2 + \frac{d-c}{(1+exp(b(x-e)))^f}}
\bold{b6:} \deqn{f(x) = c + k \cdot x + \frac{d-c}{(1+exp(b(x-e)))^f}}
\bold{b5:} \deqn{f(x) = c + \frac{d-c}{(1+exp(b(x-e)))^f}}
\bold{b4:} \deqn{f(x) = c + \frac{d-c}{1+exp(b(x-e))}}
\bold{b3:} \deqn{f(x) = \frac{d}{1+exp(b(x-e))}} 
\bold{expGrowth}: \deqn{f(x) = a \cdot exp(b \cdot x) + c}\cr

The following mechanistic models are implemented:\cr\cr
\bold{mak2}: \deqn{F_n = F_{n-1} + k \cdot log \left(1 + \left(\frac{F_{n-1}}{k}\right)\right) + Fb}
\bold{mak3 & mak3n}: \deqn{F_n = F_{n-1} + k \cdot log \left(1 + \left(\frac{F_{n-1}}{k}\right)\right) + (slope \cdot n + Fb)} 
\bold{chag}: \deqn{F_n = 1 - (1 - F_{n-1}) \cdot \left(\frac{1 - b \cdot  F_{n-1}}{1 + (a - 2) \cdot b \cdot F_{n-1}}\right)^\frac{1}{a - 1}}

\code{mak2} and \code{mak3} are two mechanistic models developed by Gregory Boggy (see references). \code{chag} was developed by Alexander Chagovetz and James Keener.

The mechanistic models are a completely different approach in that the response value (Fluorescence) is not a function of the predictor value (Cycles), but a function of the preceeding response value, that is, \eqn{F_n = f(F_{n-1})}. These are also called 'recurrence relations' or 'iterative maps'.

The implementation of these models in the 'qpcR'' package is the following:\cr
1) In case of \code{mak2} or \code{mak3}, all cycles up from the second derivative maximum of a four-parameter log-logistic model (l4) are chopped off. This is because these two models do not fit to a complete sigmoidal curve. For \code{chag} and \code{mak3n}, the sigmoidal curve is also rescaled between [0, 1].\cr 
2) A grid of sensible starting values is created for all parameters in the model.\cr
3) For each combination of starting parameters, the model is fit by nonlinear least-squares (Levenberg-Marquardt by default).\cr
4) The acquired parameters are collected in a parameter matrix together with the residual sum-of-squares (RSS) of the fit.\cr
5) The parameter combination is selected that delivered the lowest RSS.\cr
6) These parameters are transferred to \code{\link{pcrfit}}, and the data is refitted.\cr
7) Parameter \code{D0} can be used directly to calculate expression ratios, hence making the use of threshold cycles and efficiencies expendable.\cr
Some parameters for this function can be set in \code{\link{parMAK}}, see also examples there and in \code{\link{pcrfit}}.\cr
 
The functions are defined as a list containing the following items:\cr

\code{$expr}        the function as an expression for the fitting procedure.\cr
\code{$fct}         the function defined as \code{f(x, parm)}.\cr
\code{$ssfct}       the self-starter function.\cr
\code{$d1}          the first derivative function.\cr
\code{$d2}          the second derivative function.\cr
\code{$inv}         the inverse function.\cr
\code{$expr.grad}   the function as an expression for gradient calculation.\cr
\code{$inv.grad}    the inverse functions as an expression for gradient calculation.\cr
\code{$parnames}    the parameter names.\cr
\code{$name}        the function name.\cr
\code{$type}        the function type as a character string.\cr  
}

\author{
Andrej-Nikolai Spiess
}

\references{
\bold{4-parameter logistic}:\cr
Validation of a quantitative method for real time PCR kinetics.\cr
Liu W & Saint DA.\cr
\emph{Biochem Biophys Res Commun} (2002), \bold{294}:347-53. 

Standardized determination of real-time PCR efficiency from a single reaction set-up.\cr
Tichopad A, Dilger M, Schwarz G & Pfaffl MW.\cr
\emph{Nucleic Acids Res} (2003), \bold{31}:e122. 

Sigmoidal curve-fitting redefines quantitative real-time PCR with the prospective of developing automated high-throughput applications.\cr
Rutledge RG.\cr
\emph{Nucleic Acids Res} (2004), \bold{32}:e178. 

A kinetic-based sigmoidal model for the polymerase chain reaction and its application to high-capacity absolute quantitative real-time PCR.\cr
Rutledge RG & Stewart D.\cr
\emph{BMC Biotechnol} (2008), \bold{8}:47.

Evaluation of absolute quantitation by nonlinear regression in probe-based real-time PCR.\cr
Goll R, Olsen T, Cui G & Florholmen J.\cr
\emph{BMC Bioinformatics} (2006), \bold{7}:107

Comprehensive algorithm for quantitative real-time polymerase chain reaction.\cr
Zhao S & Fernald RD.\cr
\emph{J Comput Biol} (2005), \bold{12}:1047-64. 

\bold{4-parameter log-logistic; 5-parameter logistic/log-logistic}:\cr
qpcR: an R package for sigmoidal model selection in quantitative real-time polymerase chain reaction analysis.\cr
Ritz C & Spiess AN.\cr
\emph{Bioinformatics} (2008), \bold{24}:1549-51. 

Highly accurate sigmoidal fitting of real-time PCR data by introducing a parameter for asymmetry.\cr
Spiess AN, Feig C & Ritz C.\cr
\emph{BMC Bioinformatics} (2008), \bold{29}:221. 

\bold{exponential model}:\cr
Standardized determination of real-time PCR efficiency from a single reaction set-up.\cr
Tichopad A, Dilger M, Schwarz G & Pfaffl MW.\cr
\emph{Nucleic Acids Research} (2003), \bold{31}:e122.

Comprehensive algorithm for quantitative real-time polymerase chain reaction.\cr
Zhao S & Fernald RD.\cr
\emph{J Comput Biol} (2005), \bold{12}:1047-64.

\bold{mak2, mak3, mak3n}:\cr
A Mechanistic Model of PCR for Accurate Quantification of Quantitative PCR Data.\cr
Boggy GJ & Woolf PJ.\cr
\emph{PLoS ONE} (2010), \bold{5(8)}: e12355.

\bold{chag}:\cr
Kinetic models of qPCR as applied to the problem of low copy numbers.\cr
Chagovetz AM & Keener JP.\cr
Poster presented at the qPCR 2011 Symposium in Munich.\cr
Can be downloaded at \url{http://www.dr-spiess.de/qpcR/other/Chago.pdf}.
}

\examples{
m1 <- pcrfit(reps, 1, 2, b3)
m2 <- pcrfit(reps, 1, 2, b5)
m3 <- pcrfit(reps, 1, 2, l6)
m4 <- pcrfit(reps, 1, 2, l7)

## get the second derivative
## curve of m2
d2 <- b5$d2(m2$DATA[, 1], coef(m2))
plot(m2)
lines(d2, col = 2)  
}

\keyword{models}
\keyword{nonlinear}
