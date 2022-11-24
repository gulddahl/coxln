#' Create a covariance function

#' @description
#' Creates an isotropic covariance (or, more precisely, correlation)
#' function for use on linear networks.

#' @details
#' The functions returns a covariance (or correlation) function.
#' If par=NULL, this is a function of (par,r) where par is a
#' parameter vector, while if par is a vector it returns a function
#' of r. A supplied parameter vector should have the correct number
#' of entries for the chosen covariance function. The currently
#' implemented covariance functions are:
#'
#' * Exponential covariance function: c0(r) = exp(-sr),
#' with a scale parameter s>0. Set type="expcov" and par=c(s).
#' * Covariance function with gamma Bernstein CDF:
#' c0(r) = (1+r/phi)^(-tau), with shape parameter tau and rate
#' parameter phi. Set type="gamma" and par = c(tau,phi).
#' * Covariance function with inverse gamma Bernstein CDF:
#' c0(r) = 2(r phi)^(tau/2)K_tau(2 sqrt(r phi)/(Gamma(tau)),
#' with shape parameter tau and scale
#' parameter phi. Set type="invgamma" and par = c(tau,phi).
#' * Covariance function with generalized inverse Gaussian Bernstein
#' CDF: c0(r) = (1+2r/psi)^(-lambda/2)(K_lambda(sqrt((2r+psi)chi))/
#' (K_lambda(sqrt(psi chi))
#' with parameters (chi,psi,lambda); set type="gig" and
#' par = c(chi,psi,lambda).

#' @param type The types implemented are "expcov", "gamma", "invgamma"
#' or "gig"; see the details below.
#' @param par if par=NULL (the default behavior), the parameters
#' can be specified later; otherwise this should be a vector containing
#' the correct number of parameters; see the details below.
#' @returns A function of distance r, or if par=NULL, a function of
#' (par,r) where par is the parameter vector.

#' @examples
#' # Exponential covariance function with scale parameter 1
#' c0 = covfunctypes("expcov",1)
#' curve(c0(x),from=0,to=5)
#'
#' # Exponential covariance function with unspecified parameter
#' c0 = covfunctypes("expcov")
#' curve(c0(1,x),from=0,to=5)  # scale=1 specified here
#'
#' # Covariance function with gamma Bernstein CDF
#' c0 = covfunctypes("gamma",par=c(1,0.5))  # tau=1 and phi=0.5
#' curve(c0(x),from=0,to=5)
#'
#' # Covariance function with inverse gamma Bernstein CDF
#' c0 = covfunctypes("invgamma",par=c(1,2))  # tau=1 and phi=2
#' curve(c0(x),from=0,to=5)
#'
#' # Covariance function with generalized inverse Gaussian Bernstein CDF
#' c0 = covfunctypes("gig",par=c(1,1,1))  # (chi,psi,lambda) = (1,1,1)
#' curve(c0(x),from=0,to=5)

#' @export

covfunctypes = function(type,par=NULL){
  if (type=="expcov") covfunc = function(par,r) exp(-par[1]*r)
  else if (type=="gamma") covfunc = function(par,r) (1+r/par[2])^(-par[1])
  else if (type=="invgamma") covfunc = function(par,r) ifelse(r>0,2*par[2]^par[1]/gamma(par[1])*(r/par[2])^(par[1]/2)*besselK(2*sqrt(r*par[2]),par[1]),1)
  else if (type=="gig") covfunc = function(par,r) (1+2*r/par[2])^(-par[3]/2)*besselK(sqrt((2*r+par[2])*par[1]),par[3])/besselK(sqrt(par[2]*par[1]),par[3])
  else stop("unknown type of covariance function")
  if (!is.null(par)) {
    covfunc2 = function(r) covfunc(par,r)
    return(covfunc2)
  } else return(covfunc)
}


#' Create a pair correlation function

#' @description
#' This creates a theoretical pair correlation function for a
#' LGCP, ICP or PCPP model with a chosen covariance function
#' for the underlying Gaussian process.

#' @details
#' This function creates a pair correlation function with inputs
#' par and r, where par is a vector of parameters used both in the
#' chosen point process, and the chosen covariance function, and r
#' is the distance at with the pair correlation function is evaluated.
#' The parameters in par depend on the point process and the covariance
#' function, where the first parameters are used by the choice of
#' point process as follow:
#' * "lgcp": par = (sigma, ...)
#' * "icp": par = (sigma, h, ...)
#' * "pcpp": par = (h, ...)
#'
#' The remaining parameters in the par vector is used by the
#' covariance function (see the covfunctypes function). Note that
#' the number of parameters in the par vector should fit the total
#' number of parameters used by the point process and the covariance
#' function.

#' @param covtype The type of covariance function used;
#' currently available: "expcov", "gamma", "invgamma" and "gig";
#' see details under the function covfunctypes.
#' @param transform The type of point process used; "lgcp" for
#' log Gaussian Cox process, "icp" for interrupted Cox process, and
#' "pcpp" for permanental Cox process.
#' @returns A function with inputs par and r.

#' @examples
#' # pcf for LGCP with exponential covariance function
#' pcf = paircorfunc("expcov","lgcp")
#' parameters = c(sigma=1,s=1)
#' curve(pcf(parameters,x),from=0,to=5,ylab="pcf",xlab="r")
#'
#' # pcf for ICP with gamma Bernstein CDF
#' pcf = paircorfunc("gamma","icp")
#' parameters = c(sigma=1,h=1,tau=1,phi=1)
#' curve(pcf(parameters,x),from=0,to=5,ylab="pcf",xlab="r")
#'
#' # pcf for PCPP with inverse gamma Bernstein CDF
#' pcf = paircorfunc("invgamma","pcpp")
#' parameters = c(h=1,tau=1,phi=1)
#' curve(pcf(parameters,x),from=0,to=5,ylab="pcf",xlab="r")

#' @export

paircorfunc = function(covtype,transform){
  covfunc = covfunctypes(covtype)
  if (transform=="lgcp") g = function(par,r) exp(par[1]^2*covfunc(par[-1],r))
  if (transform=="icp") g = function(par,r) ((1+par[1]^2)^2/((1+par[1]^2)^2-par[1]^(2*2)*covfunc(par[-c(1:2)],r)^2))^(par[2]/2)
  if (transform=="pcpp") g = function(par,r) 1+covfunc(par[-1],r)^2/(par[1]/2)
  return(g)
}
