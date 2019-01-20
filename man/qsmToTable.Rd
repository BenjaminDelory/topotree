\name{qsmToTable}
\alias{qsmToTable}
\title{Import QSM Files Into a Single Data Frame}
\description{Import QSM files located in a folder into a single data frame.  
}
\usage{
qsmToTable(inputqsm)
}
\arguments{
  \item{inputqsm}{A character string specifying the path to the folder containing the QSM files.}
}
\value{
Returns a qsmToTable object. Each line of the data frame gives the attributes of one segment. The data frame includes the following columns:

\bold{file}: the name of each QSM file located in \code{inputqsm}.

\bold{branch}: the branch identification number.

\bold{parentbranch}: the identification number of the parent branch.

\bold{order}: the branching order.

\bold{bran}: the branching state (\sQuote{true} denotes a branching link while \sQuote{false} denotes a continuing link).

\bold{apic}: the apical state (\sQuote{true} denotes a segment ending a branch).

\bold{x1}: the X coordinate of the node starting the segment.

\bold{y1}: the Y coordinate of the node starting the segment.

\bold{z1}: the Z coordinate of the node starting the segment.

\bold{x2}: the X coordinate of the node ending the segment.

\bold{y2}: the Y coordinate of the node ending the segment.

\bold{z2}: the Z coordinate of the node ending the segment.

\bold{length}: the length of the segment.

\bold{blength}: the cumulative length measured along each branch.

\bold{geodesic}: the geodesic distance measured along the branches.
}
\author{Benjamin M. Delory
}