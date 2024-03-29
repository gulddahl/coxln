% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distR.R
\name{distR}
\alias{distR}
\title{Resistance distance between two points in a linear network}
\usage{
distR(segu, tu, segv, tv, L, ILM = solve(CalcLaplacianMatrix(L)))
}
\arguments{
\item{segu}{The line segment number of point u.}

\item{tu}{The position of u on the line segment (a value between 0 and 1)}

\item{segv}{The line segment number of point v.}

\item{tv}{The position of v on the line segment (a value between 0 and 1)}

\item{L}{A linear network of class linnet from the spatstat package.}

\item{ILM}{For internal use.}
}
\value{
The resistance distance between u and v.
}
\description{
Calculates the resistance distance between two points u and v in a
linear network.
}
\examples{
# Distance between midpoints in segments number 1 and 2
distR(1,0.5,2,0.5,simplenet)
}
\seealso{
\link{pairdistR} for calculation of the resistance
metric on one or two point pattern on a linear network,
\link{CalcLaplacianMatrix} for calculation of the Laplacian matrix
}
