% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simGaus.R
\name{simGausExpLNRoot}
\alias{simGausExpLNRoot}
\title{Simulate Gaussian Random Field with exponential covariance function on tree-shaped linear network}
\usage{
simGausExpLNRoot(
  pos,
  s,
  mu = 0,
  sigma = 1,
  transform = "none",
  h = 1,
  rho = 1,
  orderV = makeorderV(as.linnet(pos))
)
}
\arguments{
\item{pos}{A point pattern on L defining the discretisation of the network.}

\item{s}{The scale parameter used in the exponential covariance function.}

\item{mu}{The mean of the GRF. Only used if transform=="none".}

\item{sigma}{The standard deviation used in the GRF. Ignored if transform=="pcpp".}

\item{transform}{The transform applied to the GRF. "none" gives an untransformed GRF,
"lgcp" gives a log GRF (used in a log Gaussian Cox process), "icp" gives the field
used in an interupted Cox process, and "pcpp" gives the field used in a permanental
Cox point process.}

\item{h}{The number of independent GRFs used in "icp" or "pcpp"; ignored if transform
is "none" and "lgcp".}

\item{rho}{The mean number of points per unit length of network when GRF is used
for Cox point processes; ignored if transform=="none".}

\item{orderV}{The order on the vertices used for defining the order of the
line segments used in the simulation; default is the order given by the function makeorderV.
Warning: Using other orders may result in a simulation with the wrong distribution.}
}
\value{
A simulation of a (transformed) Gaussian process with class linfun from spatstat.
}
\description{
Simulates a Gaussian Random Field (GRF) with an exponential covariance function
defined on a tree-shaped linear network using sequential simulation on the line segments
one at a time starting from the root. Also includes the possibility of
transforming the simulation for use in various cox point processes.
}
\examples{
# Gaussian process on network from dendrite data
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
X = simGausExpLNRoot(pos,s=0.01,mu=0,sigma=1,transform="none")
plot(X)

# simulation of intensity for LGCP
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
X = simGausExpLNRoot(pos,s=0.01,sigma=1,transform="lgcp",rho=1)
plot(X)

# simulation of intensity for ICP
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
X = simGausExpLNRoot(pos,s=0.01,sigma=1,transform="icp",h=1,rho=1)
plot(X)

# simulation of intensity for PCPP
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
X = simGausExpLNRoot(pos,s=0.01,transform="pcpp",h=1,rho=1)
plot(X)
}
\seealso{
\link{simCPExpLNRoot} for simulation of Cox processes using GRFs,
\link{simGausLNRoot} for simulation of GRFs with arbitrary
covariance functions, \link{simGausLNDisc} for simulation of GRFs
on arbitrary linear networks
}
