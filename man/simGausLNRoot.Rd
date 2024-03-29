% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simGaus.R
\name{simGausLNRoot}
\alias{simGausLNRoot}
\title{Simulate Gaussian Random Field with arbitrary covariance function on tree-shaped linear network}
\usage{
simGausLNRoot(
  pos,
  simalgo,
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

\item{simalgo}{A vector of values used as parameters in the underlying
exponential covariance functions; can be created conveniently by the function
simalgotypes to match a particular covariance function.}

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
Simulates a Gaussian random Field with arbitrary covariance function on a
tree-shaped linear network using average of simulations with
exponential covariance functions with parameters drawn from
a given Bernstein density.
}
\details{
This function simulates a Gaussian random Field with
arbitrary covariance function on a
tree-shaped linear network using average of simulations with
exponential covariance functions with parameters drawn from
a given Bernstein density.
Each of the simulations with exponential
covariance function are simulated with a sequential simulation
algorithm simulating one line segment at a time starting at the root
of the tree. Note that the algorithm is only approximate, but a central limit
theorem ensures it converges to the right distribution, when the number of
simulations from the Bernstein density increases.
}
\examples{
# Gaussian process on network with gamma Bernstein density
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
simalgo = simalgotypes(param=c(5,5),type="gammabd",nsim=50)
X = simGausLNRoot(pos,simalgo,mu=0,sigma=1,transform="none")
plot(X)

# simulation of intensity for LGCP with gamma Bernstein density
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
simalgo = simalgotypes(param=c(5,5),type="gammabd",nsim=50)
X = simGausLNRoot(pos,simalgo,sigma=1,transform="lgcp",rho=1)
plot(X)

# simulation of intensity for ICP with inverse gamma Bernstein density
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
simalgo = simalgotypes(param=c(5,5),type="invgammabd",nsim=50)
X = simGausLNRoot(pos,simalgo,sigma=1,transform="icp",h=1,rho=1)
plot(X)

# simulation of intensity for PCPP with generalized inverse Gaussian Bernstein density
pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
simalgo = simalgotypes(param=c(5,5,5),type="gigbd",nsim=50)
X = simGausLNRoot(pos,simalgo,transform="pcpp",h=1,rho=1)
plot(X)
}
\seealso{
\link{simCPLNRoot} for simulation of Cox processes using GRFs,
\link{simGausExpLNRoot} for simulation of GRFs with exponential
covariance functions, \link{simGausLNDisc} for simulation of GRFs
on arbitrary linear networks, \link{simalgotypes} for simulation of Bernstein
distributions, \link{covfunctypes} for covariance functions, \link{paircorfunc}
for pair correlation functions
}
