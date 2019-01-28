\name{perhomology}
\alias{perhomology}
\title{Topological Analysis using Persistent Homology}
\description{Computing the birth and death of homology group components (path-connected components) using tree architecture data found in QSM files.  
}
\usage{
perhomology(x)
}
\arguments{
  \item{x}{A qsmToTable object. Mandatory.}
}
\details{
In mathematics, homology refers to homology groups. Because tree branches are path-connected components, only zero order homology (H0) branches are considered here.

The persistence of each H0 branch is computed using a geodesic distance function. More information can be found in \cite{Li et al (2017)}.
}
\value{
Returns a matrix. The matrix is a barcode object and has the following columns: 

\bold{dimension}: the homology group (0 for path-connected components).

\bold{birth}: the birth of an H0 branch along the distance function.

\bold{death}: the death of an H0 branch along the distance function.
}
\references{
Li M., Duncan K., Topp C.N., Chitwood D.H. (2017) Persistent homology and the branching topologies of plants, \emph{American Journal of Botany}, DOI: 10.3732/ajb.1700046.
}
\author{Benjamin M. Delory
}
\seealso{\code{\link{qsmToTable}}, \code{\link{plot.barcode}}.}