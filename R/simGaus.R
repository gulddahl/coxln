#' Simulate Gaussian process on linear network

#' @description
#' Simulates a Gaussian process (or transformed Gaussian process) on a linear network with
#' an arbitrary covariance function depending on the geodesic or resistance metric by
#' discretising the linear network and using the multivariate normal distribution on the grid.

#' @param pos A point pattern on the linear network representing the grid on which the process
#' is simulated. This point pattern can be made with the function makepos.
#' @param covfunc A covariance (or rather correlation) function, i.e. a function of a single parameter r representing
#' distance. Note that the covariance function should be valid, i.e. it has to be positive
#' (semi) definite.
#' @param mu The mean of the Gaussian process. If it is a number, then a constant mean is used,
#' while if it is a vector of the same length as the number of points in pos, then different means
#' are used at each grid point. Defaults to 0, and only used if transform="none".
#' @param sigma The standard deviation of the Gaussian process. Defaults to 1; and only used
#' if transform equal "none", "lgcp" or "icp"
#' @param metric The metric used by the covariance function. Either "G" for geodesic or
#' shortest path metric, or "R" for resistance metric. Defaults to "G".
#' @param transform Type of transformation used on the Gaussian process. If transform="none"
#' a Gaussian process is returned, and if transform is equal to "lgcp", "icp" or "pcpp"
#' the random intensity of a log Gaussian Cox process, interrupted Cox process or
#' permanental Cox point process is returned.
#' @param h The number of Gaussian processes simulated. Only used if transform is equal
#' to "icp" or "pcpp". Defaults to 1.
#' @param rho The intensity parameter used if transform is equal to "lgcp", "icp", or "pcpp".
#' @returns A simulation of a (transformed) Gaussian process with class linfun from spatstat.

#' @seealso [simCPLNDisc] for simulation of Cox processes using GRFs,
#' [simGausExpLNRoot] or [simGausLNDisc] for fast simulation of GRFs
#' with exponential or arbitrary covariance functions on tree-shaped
#' linear networks, [covfunctypes] for covariance functions, [paircorfunc]
#' for pair correlation functions

#' @examples
#' # Gaussian process with exponential covariance function and geodesic metric
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,mu=0,sigma=1)
#' plot(X)
#'
#' # Gaussian process with exponential covariance function and resistance metric
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,mu=0,sigma=1,metric="R")
#' plot(X)
#'
#' # Gaussian process with covariance function with gamma Bernstein CDF
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("gammabd",c(1,1))
#' X = simGausLNDisc(pos,covfunc,mu=0,sigma=1)
#' plot(X)
#'
#' # Simulation of random intensity for log Gaussian Cox process
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,sigma=1,transform="lgcp",rho=10)
#' plot(X)
#'
#' # Simulation of random intensity for interrupted Cox process
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,sigma=1,transform="icp",h=1,rho=10)
#' plot(X)
#'
#' # Simulation of random intensity for permanental Cox point process
#' pos = makepos(simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,transform="pcpp",h=1,rho=10)
#' plot(X)

#' @export

simGausLNDisc = function(pos,covfunc,mu=0,sigma=1,metric="G",transform="none",h=1,rho=1){
  if (metric=="G") covmat = covfunc(pairdist.lpp(pos))
  else if (metric=="R") covmat = covfunc(pairdistR(pos))
  else stop("Unknown metric")
  if (min(eigen(covmat)$values)+10^-6<0) stop("Covariance function is not positive (semi) definite on L.")
  if (transform=="none") Y = sigma*as.vector(rmvnorm(1,sigma=covmat)) + mu
  else if (transform=="lgcp") Y = rho*exp(sigma*as.vector(rmvnorm(1,sigma=covmat))-sigma^2/2)
  else if (transform=="icp"|transform=="pcpp") {
    Y = rep(0,npoints(pos))
    for (i in 1:h) Y = Y + (as.vector(rmvnorm(1,sigma=covmat)))^2
    if (transform=="icp") Y = rho*exp(-sigma^2*Y)*(1+2*sigma^2)^(h/2)
    if (transform=="pcpp") Y = rho*Y/h
  }
  pos = setmarks(pos,Y)
  f = nnfun.lpp(pos,value="mark")
  return(f)
}


#' Check whether L is a connected tree

#' @description
#' Checks whether L is a connected tree.

#' @param L A linear network
#' @returns Returns TRUE is L is both connected and a tree; returns FALSE otherwise.

#' @seealso [simCPLNRoot] or [simCPExpLNRoot] for simulation of Cox processes
#' specifically on trees, [simGausLNRoot] or [simGausExpLNRoot] for simulation of
#' Gaussian Random Fields specifically on trees

#' @examples
#' # Check whether simplenet from spatstat is a connected tree
#' is.tree(simplenet)
#'
#' # Check whether the network used in the dendrite data is a connected tree
#' is.tree(as.linnet(dendrite))

#' @export

is.tree = function(L) return(is.connected.linnet(L)&npoints(L$vertices)==L$lines$n+1)


#' Make order on vertices

#' @description
#' Makes an order of the vertices in a linear network (defined by the smallest
#' number of line segments from the vertex internally numbered as 1).

#' @param L A linear network
#' @returns Returns a vector of integers representing the order of each vertex.

#' @seealso [makeorderL] for orders on line segments, [simGausExpLNRoot] or
#' [simGausLNRoot] for simulation of GRFs using orders, [simCPExpLNRoot] or
#' [simCPLNRoot] for simulation of Cox processes using orders

#' @examples
#' # Make order on the vertices from simplenet
#' makeorderV(simplenet)
#'
#' # Make order on the vertices in the network used in the dendrite data
#' makeorderV(as.linnet(dendrite))

#' @export

makeorderV = function(L){
  order = Lnorder1 = c(1,rep(0,npoints(L$vertices)-1))
  i = 1
  repeat{
    i = i + 1
    Lnorder1 = L$m%*%Lnorder1
    order[as.vector(Lnorder1>0&order==0)] = i
    if (sum(order==0)==0) break
  }
  return(order)
}


#' Make order on line segments

#' @description
#' Makes an order of the line segments in a linear network by giving the same order
#' to a line segments as the lowest order of the vertices it connects to.

#' @param L A linear network
#' @param orderV A vector of orders on the vertices; created by the functino makeorderV by
#' default.
#' @returns Returns a vector of integers representing the order of each line segment.

#' @seealso [makeorderV] for orders on vertices, [simGausExpLNRoot] or
#' [simGausLNRoot] for simulation of GRFs using orders, [simCPExpLNRoot] or
#' [simCPLNRoot] for simulation of Cox processes using orders

#' @examples
#' # Make order on the vertices from simplenet
#' makeorderL(simplenet)
#'
#' # Make order on the vertices in the network used in the dendrite data
#' makeorderL(as.linnet(dendrite))

#' @export

makeorderL = function(L,orderV=makeorderV(L)){
  orderLmat = cbind(orderV[L$from],orderV[L$to])
  apply(orderLmat,1,min)
}


#' Make simulations for simGausExpLNRoot

#' @description This is an internal function for making the simulations
#' used in the function simGausExpLNRoot.
#'
#' @param pos A point pattern on L defining the discretisation of the network.
#' @param s The scale parameter used in the exponential covariance function.
#' @param mu The mean of the GRF.
#' @param sigma The standard deviation used in the GRF.
#' @param orderL The order on the line segments used in the simulation.
#' @param orderV The order on the vertices used for defining the order of the
#' line segments used in the simulation.
#' @returns A simulation of a (transformed) Gaussian process with class linfun from spatstat.

simY = function(pos,s,mu,sigma,orderL,orderV){
  L = as.linnet(pos)
  l = lengths_psp(L$lines)
  Y = rep(0,npoints(pos))
  Yvertices = c(rnorm(1),rep(0,npoints(L$vertices)-1))
  for (i in 1:max(orderL)){
    for (j in (1:length(orderL))[orderL==i]){
      if (orderV[L$from[j]]<orderV[L$to[j]]){
        posj = pos$data$tp[pos$data$seg==j]
        simj = c(Yvertices[L$from[j]],rep(0,length(posj)-1))
        if (length(posj)<2) stop("pos should contain duplicated points at vertices (i.e. use makepos function with duplicate=TRUE)")
        for (k in 2:length(posj)){
          simj[k] = rnorm(1,exp(-l[j]*(posj[k]-posj[k-1])*s)*simj[k-1],
                          sqrt(1-exp(-2*l[j]*(posj[k]-posj[k-1])*s)))
        }
        Yvertices[L$to[j]] = simj[length(posj)]
        Y[pos$data$seg==j]=simj
      } else {
        posj = pos$data$tp[pos$data$seg==j]
        simj = c(rep(0,length(posj)-1),Yvertices[L$to[j]])
        for (k in (length(posj)-1):1){
          simj[k] = rnorm(1,exp(-l[j]*(posj[k+1]-posj[k])*s)*simj[k+1],
                          sqrt(1-exp(-2*l[j]*(posj[k+1]-posj[k])*s)))
        }
        Yvertices[L$from[j]] = simj[1]
        Y[pos$data$seg==j]=simj
      }
    }
  }
  Y = sigma*Y + mu
}


#' Simulate Gaussian Random Field with exponential covariance function on tree-shaped linear network

#' @description
#' Simulates a Gaussian Random Field (GRF) with an exponential covariance function
#' defined on a tree-shaped linear network using sequential simulation on the line segments
#' one at a time starting from the root. Also includes the possibility of
#' transforming the simulation for use in various cox point processes.

#' @param pos A point pattern on L defining the discretisation of the network.
#' @param s The scale parameter used in the exponential covariance function.
#' @param mu The mean of the GRF. Only used if transform=="none".
#' @param sigma The standard deviation used in the GRF. Ignored if transform=="pcpp".
#' @param transform The transform applied to the GRF. "none" gives an untransformed GRF,
#' "lgcp" gives a log GRF (used in a log Gaussian Cox process), "icp" gives the field
#' used in an interupted Cox process, and "pcpp" gives the field used in a permanental
#' Cox point process.
#' @param h The number of independent GRFs used in "icp" or "pcpp"; ignored if transform
#' is "none" and "lgcp".
#' @param rho The mean number of points per unit length of network when GRF is used
#' for Cox point processes; ignored if transform=="none".
#' @param orderV The order on the vertices used for defining the order of the
#' line segments used in the simulation; default is the order given by the function makeorderV.
#' Warning: Using other orders may result in a simulation with the wrong distribution.
#' @returns A simulation of a (transformed) Gaussian process with class linfun from spatstat.

#' @seealso [simCPExpLNRoot] for simulation of Cox processes using GRFs,
#' [simGausLNRoot] for simulation of GRFs with arbitrary
#' covariance functions, [simGausLNDisc] for simulation of GRFs
#' on arbitrary linear networks

#' @examples
#' # Gaussian process on network from dendrite data
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simGausExpLNRoot(pos,s=0.01,mu=0,sigma=1,transform="none")
#' plot(X)
#'
#' # simulation of intensity for LGCP
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simGausExpLNRoot(pos,s=0.01,sigma=1,transform="lgcp",rho=1)
#' plot(X)
#'
#' # simulation of intensity for ICP
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simGausExpLNRoot(pos,s=0.01,sigma=1,transform="icp",h=1,rho=1)
#' plot(X)
#'
#' # simulation of intensity for PCPP
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' X = simGausExpLNRoot(pos,s=0.01,transform="pcpp",h=1,rho=1)
#' plot(X)

#' @export

simGausExpLNRoot = function(pos,s,mu=0,sigma=1,transform="none",h=1,rho=1,orderV=makeorderV(as.linnet(pos))){
  L = as.linnet(pos)
  if(is.tree(L)==FALSE) stop("The linear network is not a connected tree.")
  orderL = makeorderL(L,orderV)
  if (transform=="none") Y = simY(pos,s,mu,sigma,orderL,orderV)
  else if (transform=="lgcp") Y = rho*exp(simY(pos,s,mu=-sigma^2/2,sigma,orderL,orderV))
  else if (transform=="icp"|transform=="pcpp") {
    Y = rep(0,npoints(pos))
    for (i in 1:h) Y = Y + simY(pos,s,0,1,orderL,orderV)^2
    if (transform=="icp") Y = rho*exp(-sigma^2*Y)*(1+2*sigma^2)^(h/2)
    if (transform=="pcpp") Y = rho*Y/h
  } else {
    stop("Invalid transform: use none, lgcp, icp or pcpp.")
  }
  pos = setmarks(pos,Y)
  f = nnfun(pos,value="mark")
  return(f)
}


#' Simulate Gaussian Random Field with arbitrary covariance function on tree-shaped linear network

#' @description
#' Simulates a Gaussian random Field with arbitrary covariance function on a
#' tree-shaped linear network using average of simulations with
#' exponential covariance functions with parameters drawn from
#' a given Bernstein density.

#' @details
#' This function simulates a Gaussian random Field with
#' arbitrary covariance function on a
#' tree-shaped linear network using average of simulations with
#' exponential covariance functions with parameters drawn from
#' a given Bernstein density.
#' Each of the simulations with exponential
#' covariance function are simulated with a sequential simulation
#' algorithm simulating one line segment at a time starting at the root
#' of the tree. Note that the algorithm is only approximate, but a central limit
#' theorem ensures it converges to the right distribution, when the number of
#' simulations from the Bernstein density increases.

#' @param pos A point pattern on L defining the discretisation of the network.
#' @param simalgo A vector of values used as parameters in the underlying
#' exponential covariance functions; can be created conveniently by the function
#' simalgotypes to match a particular covariance function.
#' @param mu The mean of the GRF. Only used if transform=="none".
#' @param sigma The standard deviation used in the GRF. Ignored if transform=="pcpp".
#' @param transform The transform applied to the GRF. "none" gives an untransformed GRF,
#' "lgcp" gives a log GRF (used in a log Gaussian Cox process), "icp" gives the field
#' used in an interupted Cox process, and "pcpp" gives the field used in a permanental
#' Cox point process.
#' @param h The number of independent GRFs used in "icp" or "pcpp"; ignored if transform
#' is "none" and "lgcp".
#' @param rho The mean number of points per unit length of network when GRF is used
#' for Cox point processes; ignored if transform=="none".
#' @param orderV The order on the vertices used for defining the order of the
#' line segments used in the simulation; default is the order given by the function makeorderV.
#' Warning: Using other orders may result in a simulation with the wrong distribution.
#' @returns A simulation of a (transformed) Gaussian process with class linfun from spatstat.

#' @seealso [simCPLNRoot] for simulation of Cox processes using GRFs,
#' [simGausExpLNRoot] for simulation of GRFs with exponential
#' covariance functions, [simGausLNDisc] for simulation of GRFs
#' on arbitrary linear networks, [simalgotypes] for simulation of Bernstein
#' distributions, [covfunctypes] for covariance functions, [paircorfunc]
#' for pair correlation functions

#' @examples
#' # Gaussian process on network with gamma Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5),type="gammabd",nsim=50)
#' X = simGausLNRoot(pos,simalgo,mu=0,sigma=1,transform="none")
#' plot(X)
#'
#' # simulation of intensity for LGCP with gamma Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5),type="gammabd",nsim=50)
#' X = simGausLNRoot(pos,simalgo,sigma=1,transform="lgcp",rho=1)
#' plot(X)
#'
#' # simulation of intensity for ICP with inverse gamma Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5),type="invgammabd",nsim=50)
#' X = simGausLNRoot(pos,simalgo,sigma=1,transform="icp",h=1,rho=1)
#' plot(X)
#'
#' # simulation of intensity for PCPP with generalized inverse Gaussian Bernstein density
#' pos = makepos(as.linnet(dendrite),0.5,duplicate=TRUE)
#' simalgo = simalgotypes(param=c(5,5,5),type="gigbd",nsim=50)
#' X = simGausLNRoot(pos,simalgo,transform="pcpp",h=1,rho=1)
#' plot(X)

#' @export

simGausLNRoot = function(pos,simalgo,mu=0,sigma=1,transform="none",h=1,rho=1,orderV=makeorderV(as.linnet(pos))){
  L = as.linnet(pos)
  if(is.tree(L)==FALSE) stop("The linear network is not a tree.")
  orderL = makeorderL(L,orderV)
  s = simalgo
  if (transform=="none"|transform=="lgcp"){
    Ysum = rep(0,npoints(pos))
    for (i in 1:length(s)){
      Y = simY(pos,s[i],mu=0,sigma=1,orderL,orderV)
      Ysum = Y + Ysum
    }
    if (transform=="none") Y2 = sigma*1/sqrt(length(s))*Ysum+mu
    if (transform=="lgcp") Y2 = rho*exp(sigma*1/sqrt(length(s))*Ysum-sigma^2/2)
  }
  if (transform=="icp"|transform=="pcpp"){
    Ysum2 = rep(0,npoints(pos))
    for (j in 1:h){
      Ysum = rep(0,npoints(pos))
      for (i in 1:length(s)){
        if (transform=="icp") Y = simY(pos,s[i],mu=0,sigma=sigma,orderL,orderV)
        if (transform=="pcpp") Y = simY(pos,s[i],mu=0,sigma=1,orderL,orderV)
        Ysum = Y + Ysum
      }
      Ysum2 = Ysum2 + (1/sqrt(length(s))*Ysum)^2
    }
    if (transform=="icp") Y2 = rho*exp(-Ysum2)*(1+2*sigma^2)^(h/2)
    if (transform=="pcpp") Y2 = rho*Ysum2/h
  }
  pos = setmarks(pos,Y2)
  f = nnfun(pos,value="mark")
  return(f)
}

