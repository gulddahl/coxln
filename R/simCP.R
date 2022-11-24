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

#' @examples
#' # log Gaussian Cox process with exponential covariance function and geodesic metric
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("expcov",1)
#' X = simCPLNDisc(pos,covfunc,sigma=1,transform="lgcp",rho=5)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' plot(X,add=TRUE)
#'
#' # Interrupted Cox process with gamma Bernstein CDF and resistance metric
#' pos = makepos(spatstat.data::simplenet,50)
#' covfunc = covfunctypes("gamma",c(1,1))
#' X = simCPLNDisc(pos,covfunc,sigma=1,transform="icp",h=1,rho=5)
#' plot(attr(X,"Lambda"),style="width",col="grey",main="")
#' plot(X,add=TRUE)
#'
#' # Permanental Cox point process with inverse gamma Bernstein CDF and resistance metric
#' pos = makepos(spatstat.data::simplenet,50)
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
