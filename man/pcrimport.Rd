\name{pcrimport}
\alias{pcrimport}

\title{Simple qPCR data import function (i.e. from text files or clipboard).}

\description{
 Simple wrapper function to easily import qPCR data from the clipboard (default) or tab-delimited text files.
}

\usage{
  pcrimport(file = "clipboard", sep = "\t", header = TRUE, quote = "", 
            dec = ".", colClasses = "numeric", ...)
}

\arguments{
  \item{file}{the name of the file which the data are to be read from (full path).}
  \item{sep}{the field separator character.}
  \item{header}{a logical value indicating whether the file contains the names of the variables as its first line.}
  \item{quote}{the set of quoting characters.}
  \item{dec}{the character used in the file to denote decimal points.}
  \item{colClasses}{character. A vector of classes to be assumed for the columns.} 
  \item{...}{further arguments to be passed on to \code{\link{read.table}}.}
 }

\details{
  For a more detailed description of the arguments see \code{\link{read.table}}.
}

\value{
  A data frame containing a representation of the data in the file. 

}

\author{
  Andrej-Nikolai Spiess
}

\examples{
### paste some Excel data into the clipboard ###
\dontrun{
temp <- pcrimport()
}
### from a tab-delimited text file ###
\dontrun{
temp <- pcrimport("c:\\temp\\foo.txt")
}
}

\keyword{IO}
\keyword{file}
