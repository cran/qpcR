\name{cbind.na}
\alias{cbind.na}
\alias{rbind.na}
\alias{data.frame.na}  

\title{Create data frames from objects of unequal size or combine R objects of unequal size by Rows or Columns}  

\description{
For \code{data.frame.na}, a data frame can be created from objects of unequal length (thereby filling to maximum occuring length with \code{NA}'s).
For \code{cbind} and \code{rbind}, take a sequence of vector, matrix or data frames arguments and combine by \emph{c}olumns or \emph{r}ows,
 filling to maximum occuring length with \code{NA}'s. 
}

\usage{
cbind.na(..., deparse.level = 1)
rbind.na(..., deparse.level = 1)
data.frame.na(..., row.names = NULL, check.rows = FALSE, check.names = TRUE,
              stringsAsFactors = default.stringsAsFactors())
}

\arguments{
 \item{...}{vectors, matrices or data frames. See documentation to \code{\link{cbind}}.}
 \item{deparse.level}{integer controlling the construction of labels in the case of non-matrix-like arguments (for the default method). See documentation to \code{\link{cbind}}.}  
 \item{row.names}{NULL or a single integer or character string specifying a column to be used as row names, or a character or integer vector giving the row names for the data frame.}
 \item{check.rows}{if TRUE then the rows are checked for consistency of length and names.}
 \item{check.names}{logical. If \code{TRUE} then the names of the variables in the data frame are checked to ensure that they are syntactically valid variable names and are not duplicated. If necessary they are adjusted (by \code{\link{make.names}}) so that they are.}
 \item{stringsAsFactors}{logical: should character vectors be converted to factors? The factory-fresh default is \code{TRUE}, but this can be changed by setting \code{options(stringsAsFactors = FALSE)}.}
}

\details{
For more details, see documentation to  to \code{\link{cbind}}, \code{\link{rbind}} or \code{\link{data.frame}}.
The major difference of the three functions described here is that they work like the standard generics if sizes of vectors or matrices are equal,
 but if sizes are unequal, vectors or matrices are filled with \code{NA}'s to the maximum occuring length (or row/column number) before calling the standard generics.
 This has the advantage that in case of unequal vectors, values are not replicated to length, which changes subsequent analysis such as average calculations etc.
 Also, in case of \code{\link{data.frame}}, unequal length vectors result in a stop, which is not the case with \code{\link{data.frame.na}}.
 In case of binding unequal length matrices, \code{\link{cbind}} or \code{\link{rbind}} will stop, giving an error message about unequal row or columns sizes, respectively.
 \code{cbind.na}, \code{rbind.na} and \code{data.frame.na} always work in this case.    
}

\value{
An object with created or combined variables, filled with \code{NA}'s to equal length.  
}

\note{
For \code{cbind} and \code{rbind}, if recursively used in large size loops, it can take substantially longer than the original generics due to repeated testing for unequal length.
}

\author{
Andrej-Nikolai Spiess, taken and altered from \code{methods:::cbind}, \code{methods:::rbind} and \code{data.frame}.
}

\examples{
## binding
cbind.na(1, 1:7) # the '1' (= shorter vector) is NOT recycled but filled
cbind.na(1:8, 1:7, 1:5, 1:10) # same with many vectors
rbind.na(1:8, 1:7, 1:5, 1:10) # or in rows

a <- matrix(rnorm(20), ncol = 4) # unequal size matrices
b <- matrix(rnorm(20), ncol = 5)
cbind.na(a, b) # works, in contrast to original cbind
rbind.na(a, b) # works, in contrast to original rbind

## data frame with unequal size vectors
data.frame.na(A = 1:7, B = 1:5, C = letters[1:3], 
              D = factor(c(1, 1, 2, 2)), stringsAsFactors = FALSE) 
              
## convert a list with unequal length list items
## to a data frame
z <- list(a = 1:5, b = letters[1:3], c = matrix(rnorm(20), ncol = 2))
do.call(data.frame.na, c(z, stringsAsFactors = FALSE))
}

\keyword{classes}
\keyword{methods}
