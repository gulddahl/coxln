#' Simulate Cox process on linear network

#' @description
#' Simulates a Cox process (a log Gaussian Cox process, interrupted
#' Cox process or permanental Cox point process) on a linear
#' network with an arbitrary covariance function depending on
#' the geodesic or resistance metric by discretising the linear
#' network and using the multivariate normal distribution on the grid.

#' @param pos A point pattern on the linear network representing the grid on which the underlying
#' Gaussian process is simulated. This point pattern can be made with the function makepos.
#' @param covfunc A covariance (or rather correlation) function, i.e. a function of a single
#' parameter r representing distance. Note that the covariance function should be valid, i.e.
#' it has to be positive (semi) definite.
#' @param sigma The standard deviation of the Gaussian process. Defaults to 1; and only used
#' if transform equal "lgcp" or "icp".
#' @param metric The metric used by the covariance function. Either "G" for geodesic or
#' shortest path metric, or "R" for resistance metric. Defaults to "G".
#' @param transform Type of Cox process simulated. Possible values are "lgcp", "icp" or "pcpp"
#' corresponding to a log Gaussian Cox process, interrupted Cox process or
#' permanental Cox point process.
#' @param h The number of underlying Gaussian processes simulated. Only used if transform is equal
#' to "icp" or "pcpp". Defaults to 1.
#' @param rho The intensity of the Cox process.
#' @param savelambda If TRUE the simulation of the underlying Gaussian process is saved
#' as an attribute. Default is TRUE.
#' @returns A point process on a linear network with class lpp.

#' @seealso [simGausLNDisc] for simulation of the underlying GRF,
#' [simCPExpLNRoot] or [simCPLNRoot] for fast simulation of Cox
#' processes with exponential or arbitrary
#' covariance functions on tree-shaped linear networks,
#' [covfunctypes] for covariance functions, [paircorfunc]
#' for pair correlation functions

#' @examples
#' # log Gaussian Cox process with exponential covariance function and geodesic metric
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simCPLNDisc(pos,covfunc,sigma=1,transform="lgcp",rho=5)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' plot(X,add=TRUE)
#'
#' # Interrupted Cox process with gamma Bernstein CDF and resistance metric
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("gamma",c(1,1))
#' X = simCPLNDisc(pos,covfunc,sigma=1,transform="icp",h=1,rho=5)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' plot(X,add=TRUE)
#'
#' # Permanental Cox point process with inverse gamma Bernstein CDF and resistance metric
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("invgamma",c(1,1))
#' X = simCPLNDisc(pos,covfunc,transform="pcpp",h=1,rho=5)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' plot(X,add=TRUE)

#' @export

simCPLNDisc = function(pos,covfunc,sigma=1,metric="G",transform="lgcp",h=1,rho=1,savelambda=TRUE){
  if (transform!="lgcp"&transform!="icp"&transform!="pcpp") stop("transform must be lgcp, icp or pcpp")
  GP = simGausLNDisc(pos,covfunc,mu=0,sigma,metric,transform,h,rho)
  X = rpoislpp(GP)
  if (savelambda==TRUE) attr(X,"Lambda") = GP
  return(X)
}


#' Simulate Cox process on tree-shaped linear network

#' @description
#' Simulates a Cox process (a log Gaussian Cox process, interrupted
#' Cox process or permanental Cox point process) on a tree-shaped linear
#' network using sequential simulation on the line segments
#' one at a time starting from the root of the tree. The underlying Gaussian
#' random field is equipped with an exponential covariance function.

#' @param pos A point pattern on the linear network representing the grid on which the underlying
#' Gaussian process is simulated. This point pattern can be made with the function makepos.
#' @param s The scale parameter used in the exponential covariance function.
#' @param sigma The standard deviation of the Gaussian process. Defaults to 1 and only used
#' if transform equal "lgcp" or "icp".
#' @param transform Type of Cox process simulated. Possible values are "lgcp", "icp" or "pcpp"
#' corresponding to a log Gaussian Cox process, interrupted Cox process or
#' permanental Cox point process.
#' @param h The number of underlying Gaussian processes simulated. Only used if transform is equal
#' to "icp" or "pcpp". Defaults to 1.
#' @param rho The intensity of the Cox process.
#' @param savelambda If TRUE the simulation of the underlying Gaussian process is saved
#' as an attribute. Default is TRUE.
#' @param orderV The order on the vertices used for defining the order of the
#' line segments used in the simulation; default is the order given by the function makeorderV.
#' Warning: Using other orders may result in a simulation with the wrong distribution.
#' @returns A point process on a linear network with class lpp.

#' @seealso [simGausExpLNRoot] for simulation of the underlying GRF,
#' [simCPLNRoot] for simulation of Cox processes with arbitrary
#' covariance functions, [simCPLNDisc] for simulation of Cox processes
#' on arbitrary linear networks

#' @examples
#' # log Gaussian Cox process
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simCPExpLNRoot(pos,s=0.01,sigma=1,transform="lgcp",rho=0.05,savelambda=TRUE)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' points(X,cex=0.5)
#'
#' # Interrupted Cox process
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simCPExpLNRoot(pos,s=0.01,sigma=10,transform="icp",h=1,rho=0.05,savelambda=TRUE)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' points(X,cex=0.5)
#'
#' # Permanental Cox point process
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simCPExpLNRoot(pos,s=0.01,transform="pcpp",h=1,rho=0.05,savelambda=TRUE)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' points(X,cex=0.5)

#' @export

simCPExpLNRoot = function(pos,s,sigma=1,transform="lgcp",h=1,rho=1,orderV=makeorderV(as.linnet(pos)),savelambda=TRUE){
  if (transform!="lgcp"&transform!="icp"&transform!="pcpp") stop("transform must be lgcp, icp or pcpp")
  GP = simGausExpLNRoot(pos,s,mu=0,sigma,transform,h,rho,orderV)
  X = rpoislpp(GP)
  if (savelambda==TRUE) attr(X,"Lambda") = GP
  return(X)
}


#' Simulate Cox process with arbitrary covariance function on tree-shaped linear network

#' @description
#' Simulates a Gaussian random Field with arbitrary covariance function on a
#' tree-shaped linear network using average of simulations with
#' exponential covariance functions with parameters drawn from
#' a given Bernstein density.

#' @description
#' Simulates a Cox process (a log Gaussian Cox process, interrupted
#' Cox process or permanental Cox point process) on a tree-shaped linear
#' network by the taking the average of a number of simulations each
#' using sequential simulation on the line segments
#' one at a time starting from the root of the tree. The underlying Gaussian
#' random field is equipped with an arbitrary covariance function.

#' @details
#' This function simulates a Cox process by simulating a number of
#' Gaussian random Fields each with exponential covariance function on a
#' tree-shaped linear network with parameters drawn from
#' a given Bernstein density. Using the average of these simulations an approximate
#' simulation of a Gaussian Random Field with the covariance function corresponding
#' to the used Bernstein density is achieved for a sufficiently high number of
#' simulations, since a central limit
#' theorem ensures it converges to the right distribution, when the number of
#' simulations from the Bernstein distribution increases. The Gaussian random fields
#' are transformed using a transformation depending which type of Cox process is wanted
#' (i.e. a log Gaussian Cox process, interrupted
#' Cox process or permanental Cox point process), and the resulting field is
#' used as an intensity in a Poisson process on the linear network.

#' @param pos A point pattern on L defining the discretisation of the network.
#' @param simalgo A vector of values used as parameters in the underlying
#' exponential covariance functions; can be created conveniently by the function
#' simalgotypes to match a particular covariance function. The length of this vector
#' gives the number of underlying simulations.
#' @param sigma The standard deviation used in the GRF. Ignored if transform=="pcpp".
#' @param transform The transform applied to the GRF. "lgcp" gives a log Gaussian
#' Cox process), "icp" gives an interupted Cox process, and "pcpp" gives a permanental
#' Cox point process.
#' @param h The number of independent GRFs used in "icp" or "pcpp"; ignored if transform
#' is "lgcp".
#' @param rho The mean number of points in the Cox process per unit length of network;
#' ignored if transform=="none".
#' @param orderV The order on the vertices used for defining the order of the
#' line segments used in the simulation; default is the order given by the function makeorderV.
#' Warning: Using other orders may result in a simulation with the wrong distribution.
#' @param savelambda If TRUE the simulation of the underlying Gaussian process is saved
#' as an attribute. Default is TRUE.
#' @returns A point process on a linear network with class lpp.

#' @seealso [simGausLNRoot] for simulation of the underlying GRF,
#' [simCPExpLNRoot] for simulation of Cox processes with exponential
#' covariance functions, [simCPLNDisc] for simulation of Cox processes
#' on arbitrary linear networks, [simalgotypes] for simulation of Bernstein
#' distributions, [covfunctypes] for covariance functions, [paircorfunc]
#' for pair correlation functions

#' @examples
#' # simulation of LGCP with gamma Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5),type="gamma",nsim=50)
#' X = simCPLNRoot(pos,simalgo,sigma=1,transform="lgcp",rho=0.1)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' points(X,cex=0.5)
#'
#' # simulation of ICP with inverse gamma Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5),type="invgamma",nsim=50)
#' X = simCPLNRoot(pos,simalgo,sigma=1,transform="icp",h=1,rho=0.1)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' points(X,cex=0.5)
#'
#' # simulation of PCPP with generalized inverse Gaussian Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5,5),type="gig",nsim=50)
#' X = simCPLNRoot(pos,simalgo,transform="pcpp",h=1,rho=0.1)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' points(X,cex=0.5)

#' @export

simCPLNRoot = function(pos,simalgo,sigma=1,transform="lgcp",h=1,rho=1,orderV=makeorderV(as.linnet(pos)),savelambda=TRUE){
  if (transform!="lgcp"&transform!="icp"&transform!="pcpp") stop("transform must be lgcp, icp or ppp")
  GP = simGausLNRoot(pos,simalgo,mu=0,sigma,transform,h,rho,orderV)
  X = rpoislpp(GP)
  if (savelambda==TRUE) attr(X,"Lambda") = GP
  return(X)
}
