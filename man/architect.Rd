\name{architect}
\alias{architect}
\title{Computing Traits Describing The Global Tree Architecture}
\description{Calculation of common architectural traits from a QSM output file.}
\usage{
architect(inputqsm)}
\arguments{
  \item{inputqsm}{A qsmToTable object.}
}
\value{
Returns a data frame summarizing all architectural traits calculated for the QSM file located in /\code{inputqsm}. 
The data frame contains the following columns:

\bold{FileName}: the name of the QSM file.

\bold{TBL}: the total length of branches.

\bold{L0B}: the length of the zero-order branch.

\bold{TN0B}: the total number of zero-order branches.

\bold{TNLB}: the total number of lateral branches.

\bold{TLBL}: the total length of lateral branches.

\bold{NxLB}: the total number of lateral branches with a branching order equal to x.

\bold{LxLB}: the total length of lateral branches with a branching order equal to x.

\bold{MLxLB}: the mean length of lateral branches with a branching order equal to x.

\bold{D1LB}: the density of first-order branches on the zero-order branch.

\bold{Height}: the tree height.

\bold{Width}: the tree width.

\bold{Convexhull3D}: the volume of the 3D convex hull.

\bold{Volume}: the tree volume.

\bold{Surface}: the tree surface area.

\bold{DBH}: diameter at breast height.

}
\author{Benjamin M. Delory
}
\seealso{\code{\link{qsmToTable}}.}