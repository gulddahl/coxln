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

#' @examples
#' # Gaussian process with exponential covariance function and geodesic metric
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,mu=0,sigma=1)
#' plot(X)
#'
#' # Gaussian process with exponential covariance function and resistance metric
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,mu=0,sigma=1,metric="R")
#' plot(X)
#'
#' # Gaussian process with covariance function with gamma Bernstein CDF
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("gamma",c(1,1))
#' X = simGausLNDisc(pos,covfunc,mu=0,sigma=1)
#' plot(X)
#'
#' # Simulation of random intensity for log Gaussian Cox process
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,sigma=1,transform="lgcp",rho=10)
#' plot(X)
#'
#' # Simulation of random intensity for interrupted Cox process
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,sigma=1,transform="icp",h=1,rho=10)
#' plot(X)
#'
#' # Simulation of random intensity for permanental Cox point process
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simGausLNDisc(pos,covfunc,transform="pcpp",h=1,rho=10)
#' plot(X)

#' @export

simGausLNDisc = function(pos,covfunc,mu=0,sigma=1,metric="G",transform="none",h=1,rho=1){
  posmark = pos
  if (metric=="G") covmat = covfunc(pairdist.lpp(pos))
  else if (metric=="R") covmat = covfunc(pairdistR(pos))
  else stop("Unknown metric")
  if (min(eigen(covmat)$values)+10^-6<0) stop("Covariance function is not positive (semi) definite on L.")
  if (transform=="none") Y = sigma*as.vector(rmvnorm(1,sigma=covmat)) + mu
  else if (transform=="lgcp") Y = rho*exp(sigma*as.vector(rmvnorm(1,sigma=covmat))-sigma^2/2)
  else if (transform=="icp"|transform=="pcpp") {
    Y = rep(0,npoints(pos))
    for (i in 1:h) Y = Y + (sigma*as.vector(rmvnorm(1,sigma=covmat)))^2
    if (transform=="icp") Y = rho*exp(-Y)*(1+2*sigma^2)^(h/2)
    if (transform=="pcpp") Y = rho*Y/h
  }
  #marks(pos) = Y
  pos = spatstat.geom::setmarks(pos,Y)
  f = nnfun.lpp(pos,value="mark")
  return(f)
}
