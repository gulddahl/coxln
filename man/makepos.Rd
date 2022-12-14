% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makepos.R
\name{makepos}
\alias{makepos}
\title{Make a regular grid of points on a linear network}
\usage{
makepos(L, ppul, duplicate = FALSE)
}
\arguments{
\item{L}{A linear network of class linnet from the spatstat package.}

\item{ppul}{Number of points per unit length on each line segment.}

\item{duplicate}{If duplicate is TRUE each vertex will contain
duplicate points (one for each line segment extending
from the vertex); if duplicate is FALSE, there is only one
point on each vertex.}
}
\value{
A point pattern of class lpp from the spatstat package.
}
\description{
makepos() makes a grid of points on a linear network
by placing points on the vertices and points regularly
spaced on each line segment.
}
\examples{
makepos(spatstat.data::simplenet,10)
}
