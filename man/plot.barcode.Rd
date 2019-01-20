\name{plot.barcode}
\alias{plot.barcode}
\title{Plot the Persistence Barcode}
\description{Plot the persistence barcode of the topology of a tree.}
\usage{
\method{plot}{barcode}(x, xlab=NULL, ylab=NULL, ...)
}
\arguments{
  \item{x}{A barcode object. Mandatory.}
  \item{xlab}{A character string giving a label for the X axis.}
  \item{ylab}{A character string giving a label for the Y axis.}
  \item{...}{Additional graphical parameters.}
}
\value{
Returns a plot.
}
\author{Benjamin M. Delory
}
\seealso{\code{\link{perhomology}}, \code{\link{par}}.
}