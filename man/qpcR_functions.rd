\name{qpcR_functions}
\alias{l6}
\alias{l5}
\alias{l4}
\alias{l3}
\alias{b6}
\alias{b5}
\alias{b4}
\alias{b3}
\alias{w4}
\alias{w3}
\alias{expGrowth}
\alias{baro5}

\title{The nonlinear models implemented in qpcR}

\description{
A summary of all available models implemented in this package.
}

\usage{
l6
l5
l4
l3
b6
b5
b4
b3
w4
w3
baro5
expGrowth 
}

\details{
The following nonlinear models are implemented:\cr\cr
\bold{l6:} \deqn{f(x) = c + k \cdot log(x) + \frac{d-c}{(1+exp(b(log(x)-log(e))))^f}}
\bold{l5:} \deqn{f(x) = c + \frac{d-c}{(1+exp(b(log(x)-log(e))))^f}}
\bold{l4:} \deqn{f(x) = c + \frac{d-c}{1+exp(b(log(x)-log(e)))}}
\bold{l3:} \deqn{f(x) = \frac{d}{1+exp(b(log(x)-log(e)))}} 
\bold{b6:} \deqn{f(x) = c + k \cdot x + \frac{d-c}{(1+exp(b(log(x)-log(e))))^f}}
\bold{b5:} \deqn{f(x) = c + \frac{d-c}{(1+exp(b(x-e)))^f}}
\bold{b4:} \deqn{f(x) = c + \frac{d-c}{1+exp(b(x-e))}}
\bold{b3:} \deqn{f(x) = \frac{d}{1+exp(b(x-e))}}
\bold{w4:} \deqn{f(x) = c + (d-c) exp(-exp(b(log(x)-log(e))))}
\bold{w3:} \deqn{f(x) = d exp(-exp(b(log(x)-log(e))))}
\bold{expGrowth}: \deqn{f(x) = a * exp(b * x) + c}
\bold{baro5}: \deqn{f(x) = c + \frac{d-c}{1+fexp(b1(log(x)-log(e))) + (1-f)exp(b2(log(x)-log(e)))}} with 
              \deqn{f = \frac{1}{1 + exp((2b1b2/|b1+b2|)(log(x)-log(e)))}} 

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

\examples{
m1 <- pcrfit(reps, 1, 2, b3)
m2 <- pcrfit(reps, 1, 2, b5)
m3 <- pcrfit(reps, 1, 2, w4)

## get the second derivative
## curve of m2
d2 <- b5$d2(m2$DATA[, 1], coef(m2))
plot(m2)
lines(d2, col = 2)  
}

\keyword{models}
\keyword{nonlinear}
