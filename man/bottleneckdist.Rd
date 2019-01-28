\name{bottleneckdist}
\alias{bottleneckdist}
\title{Computing a Pairwise Bottleneck Distance Matrix}
\description{Computing a pairwise distance matrix using a bottleneck distance method for the comparison of persistence barcodes.  
}
\usage{
bottleneckdist(x, show.progress=FALSE)
}
\arguments{
  \item{x}{A list of persistence barcodes. Mandatory.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
The bottleneck distance between two persistence barcodes is calculated with the \code{bottleneck} function of the TDA package (\cite{Fasy et al (2017)}). An error message will be returned if TDA is not available.

\code{show.progress=TRUE} can only be used if \code{x} contains more than one element.
}
\value{
Returns a pairwise distance matrix.
}
\references{
Fasy B.T., Kim J., Lecci F., Maria C., Millman D.L., Rouvreau V. (2018). TDA: Statistical Tools for Topological Data Analysis. \url{https://CRAN.R-project.org/package=TDA}.
}
\author{Benjamin M. Delory
}
\seealso{\code{\link{perhomology}}.}