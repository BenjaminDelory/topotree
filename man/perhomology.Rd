\name{perhomology}
\alias{perhomology}
\title{Topological Analysis using Persistent Homology}
\description{Computing the birth and death of homology group components (path-connected components) using tree architecture data found in QSM files.  
}
\usage{
perhomology(x, show.progress=FALSE)
}
\arguments{
  \item{x}{A qsmToTable object. Mandatory.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
In mathematics, homology refers to homology groups. Because tree branches are path-connected components, only zero order homology (H0) branches are considered here.

The persistence of each H0 branch is computed using a geodesic distance function. More information can be found in \cite{Li et al (2017)}.

\code{show.progress=TRUE} can only be used if more than one root system is stored in \code{x}.
}
\value{
Returns a list of matrices. The list is a perhomology object. Each matrix of the list is a barcode object. Each element of the list is named as its corresponding QSM file. Each matrix is associated with one tree only and possesses the following columns: 

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