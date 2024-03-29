% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linearFGJ.R
\name{linearG}
\alias{linearG}
\title{Estimate G function on linear network}
\usage{
linearG(X, r, metric = "G")
}
\arguments{
\item{X}{A point pattern on a linear network.}

\item{r}{A vector of distances on which the estimate of G(r) is calculated.}

\item{metric}{The metric used to measure distances. Should be either "G"
for geodesic distance, or "R" for resistance distance.}
}
\value{
The G function with class fv from spatstat.
}
\description{
Estimates the G function for a point pattern observed on a linear network.
}
\details{
This function takes a point pattern on a linear network and estimates
the G function for a set of r-values. Both the geodesic and the
resistance metric can be used for distances
}
\examples{
# Estimate G function for spiders data
X = spiders
r = seq(0,250,length.out=50)
G = linearG(X,r,metric="G")
plot(G)

# Estimate G function for dendrite data
X = dendrite
r = seq(0,5,length.out=50)
G = linearG(X,r,metric="G")
plot(G)
}
\seealso{
\link{linearF}, \link{linearJ} and \link{linearpcfR} for other summary statistics,
\link{minussampling} for minus sampling
}
