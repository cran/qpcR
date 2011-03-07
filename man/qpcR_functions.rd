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

\title{The nonlinear models implemented in qpcR}

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
}

\details{
The following nonlinear models are implemented:\cr\cr
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

\bold{mak2} and \bold{mak3}: These are two mechanistic models developed by Gregory Boggy (see references). This is a completely different approach in that
 the response value (fluorescence) is not a function of the predictor value (cycles), but a function of the preceeding response value, that is, \eqn{F_n = f(F_{n-1})}. 
 The implementation of this model in the \emph{qpcR} package is the following:\cr
 1) Either i) a sigmoidal (l4) model is fit to the overall data or ii) the data is smoothed (depending on \code{makXpar$SS.deriv}).\cr
 2) The second derivative maximum is calculated from either i) the second derivative maximum of a 4-par sigmoidal fit or ii) from the smoothed data using two-times \code{\link{diff}}.\cr
 3) All response values above the threshold cycle are removed.\cr
 4) An exponential growth model is fit to obtain sensible starting values for \code{D0} (initial template fluorescence) and \code{Fb} (background fluorescence).\cr
 5) A grid of starting values is created from 0.0001 - 10 * \code{D0}, 0.1 - 3 in steps of 0.3 for parameter \code{k}, and \code{Fb}. These are 60 combinations in total.\cr
 6) For each combination one of the following two models are fit by nonlinear least-squares (Levenberg-Marquardt by default):
 \deqn{F_n = F_{n-1} + k \cdot log(1 + (\frac{F_{n-1}}{k})) + Fb \qquad(\mathbf{mak2})} or
 \deqn{F_n = F_{n-1} + k \cdot log(1 + (\frac{F_{n-1}}{k})) + (slope \cdot n + Fb) \qquad (\mathbf{mak3})}
 7) For each combination of starting values, the optimized parameters are collected in a parameter matrix together with the residual sum-of-squares (RSS) of the fit.\cr
 8) The combination is selected that delivered the lowest RSS.\cr
 9) These values are transferred to \code{\link{pcrfit}}, and the data is refitted with the selected model using the best parameter set from 8).\cr
 10) Parameter \code{D0} can be used directly to calculate expression ratios, hence making the use of threshold cycles and efficiencies expendable.\cr
 Some parameters for this function can be set in \code{\link{makXpar}}, see also examples there and in \code{\link{pcrfit}}.\cr
 
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
A Mechanistic Model of PCR for Accurate Quantification of Quantitative PCR Data.\cr
Boggy GJ and Woolf PJ.\cr
\emph{PLoS ONE}, \bold{5(8)}: e12355.
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
