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
#' parameter phi. Set type="gammabd" and par = c(tau,phi).
#' * Covariance function with inverse gamma Bernstein CDF:
#' c0(r) = 2(r phi)^(tau/2)K_tau(2 sqrt(r phi)/(Gamma(tau)),
#' with shape parameter tau and scale
#' parameter phi. Set type="invgammabd" and par = c(tau,phi).
#' * Covariance function with generalized inverse Gaussian Bernstein
#' CDF: c0(r) = (1+2r/psi)^(-lambda/2)(K_lambda(sqrt((2r+psi)chi))/
#' (K_lambda(sqrt(psi chi))
#' with parameters (chi,psi,lambda); set type="gigbd" and
#' par = c(chi,psi,lambda).

#' @param type The types implemented are "expcov", "gammabd", "invgammabd"
#' or "gigbd"; see the details below.
#' @param par if par=NULL (the default behavior), the parameters
#' can be specified later; otherwise this should be a vector containing
#' the correct number of parameters; see the details below.
#' @returns A function of distance r, or if par=NULL, a function of
#' (par,r) where par is the parameter vector.

#' @seealso [paircorfunc] for pair correlation functions, [simalgotypes]
#' for simulation of Bernstein distributions corresponding to covariance
#' functions, [simGausLNDisc] for simulation of GRFs on linear networks
#' with a given covariance function, [parnamespcf] for names of parameters
#' in the covariance function (or pair correlation function)

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
#' c0 = covfunctypes("gammabd",par=c(1,0.5))  # tau=1 and phi=0.5
#' curve(c0(x),from=0,to=5)
#'
#' # Covariance function with inverse gamma Bernstein CDF
#' c0 = covfunctypes("invgammabd",par=c(1,2))  # tau=1 and phi=2
#' curve(c0(x),from=0,to=5)
#'
#' # Covariance function with generalized inverse Gaussian Bernstein CDF
#' c0 = covfunctypes("gigbd",par=c(1,1,1))  # (chi,psi,lambda) = (1,1,1)
#' curve(c0(x),from=0,to=5)

#' @export

covfunctypes = function(type,par=NULL){
  if (type=="expcov") covfunc = function(par,r) exp(-par[1]*r)
  else if (type=="gammabd") covfunc = function(par,r) (1+r/par[2])^(-par[1])
  else if (type=="invgammabd") covfunc = function(par,r) ifelse(r>0,2*par[2]^par[1]/gamma(par[1])*(r/par[2])^(par[1]/2)*besselK(2*sqrt(r*par[2]),par[1]),1)
  else if (type=="gigbd") covfunc = function(par,r) (1+2*r/par[2])^(-par[3]/2)*besselK(sqrt((2*r+par[2])*par[1]),par[3])/besselK(sqrt(par[2]*par[1]),par[3])
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
#' currently available: "expcov", "gammabd", "invgammabd" and "gigbd";
#' see details under the function covfunctypes.
#' @param transform The type of point process used; "lgcp" for
#' log Gaussian Cox process, "icp" for interrupted Cox process, and
#' "pcpp" for permanental Cox process.
#' @returns A function with inputs par and r.

#' @seealso [covfunctypes] for covariance functions, [simGausLNDisc]
#' for simulation of GRFs on linear networks
#' with a given covariance function, [parnamespcf] for names of parameters
#' in the pair correlation function

#' @examples
#' # pcf for LGCP with exponential covariance function
#' pcf = paircorfunc("expcov","lgcp")
#' parameters = c(sigma=1,s=1)
#' curve(pcf(parameters,x),from=0,to=5,ylab="pcf",xlab="r")
#'
#' # pcf for ICP with gamma Bernstein CDF
#' pcf = paircorfunc("gammabd","icp")
#' parameters = c(sigma=1,h=1,tau=1,phi=1)
#' curve(pcf(parameters,x),from=0,to=5,ylab="pcf",xlab="r")
#'
#' # pcf for PCPP with inverse gamma Bernstein CDF
#' pcf = paircorfunc("invgammabd","pcpp")
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


#' Simulate a Bernstein density

#' @description
#' Makes a simulation from the Bernstein density corresponding to
#' a covariance function.

#' @details
#' This function makes a number of simulations from the Bernstein density
#' corresponding to a covariance function (see the function covfunc for details
#' on the corresponding covariance functions). This vector can be used in the
#' function simGausLNRoot for making approximate simulations of a Gaussian random
#' field with the corresponding covariance function. The currently implemented
#' types are:
#' * Gamma Bernstein density: gamma distribution with shape parameter tau and rate
#' parameter phi. Set type="gammabd" and param = c(tau,phi).
#' * Inverse Gamma Bernstein density: Inverse gamma distribution with shape
#' parameter tau and scale parameter phi. Set type="invgammabd" and param = c(tau,phi).
#' * Generalized inverse Gaussian Bernstein density: Generalized inverse Gaussian
#' distribution with with parameters (chi,psi,lambda). Set type="gigbd" and
#' param = c(chi,psi,lambda).

#' @param param Parameter vector used in the Bernstein distribution.
#' @param type The type of Bernstein distribution. Can be "gammabd",
#' "invgammabd", and "gigbd". See details.
#' @param nsim The number of simulations.
#' @returns A vector of length nsim with the simulated values.

#' @seealso [covfunctypes] for covariance functions corresponding to
#' Bernstein distributions, [simGausExpLNRoot] or
#' [simGausLNRoot] for simulation of GRFs on tree-shaped linear networks
#' with a given Bernstein distribution

#' @examples
#' # simulation of gamma Bernstein distribution
#' simalgotypes(param=c(1,1),type="gammabd",nsim=10)
#'
#' # simulation of inverse gamma Bernstein distribution
#' simalgotypes(param=c(1,1),type="invgammabd",nsim=10)
#'
#' # simulation of generalized inverse Gaussian Bernstein distribution
#' simalgotypes(param=c(1,1,1),type="gigbd",nsim=10)
#'

#' @export

simalgotypes = function(param,type,nsim){
  if (type=="gammabd") simalgo = rgamma(nsim,param[1],rate=param[2])
  else if (type=="invgammabd") simalgo = 1/rgamma(nsim,param[1],rate=param[2])
  else if (type=="gigbd") simalgo = GeneralizedHyperbolic::rgig(nsim,param=param)
  else stop("unknown type of covariance function")
  return(simalgo)
}


#' Output parameter vector of a pair correlation function

#' @description
#' Since the pair correlation function both gets parameters from the covariance function
#' and the class of point process used, the number and order of parameters in its
#' parameter vector can be confusing. This function outputs the names of parameters in the
#' parameter vector of the pair correlation function.

#' @param covtype The type of covariance function used;
#' currently available: "expcov", "gammabd", "invgammabd" and "gigbd";
#' see details under the function covfunctypes.
#' @param transform The type of point process used; "lgcp" for
#' log Gaussian Cox process, "icp" for interrupted Cox process, and
#' "pcpp" for permanental Cox process. If transform = "none" only the
#' parameters in the covariance function is returned.
#' @returns A vector of strings with the names of the parameters.

#' @seealso [paircorfunc] for pair correlation functions, [covfunctypes]
#' for covariance functions.

#' @examples
#' # parameters in pcf for LGCP with exponential covariance function
#' parnamespcf("expcov","lgcp")
#'
#' # parameters in pcf for ICP with gamma covariance function
#' parnamespcf("gammabd","icp")
#'
#' # parameters in covariance function with generalized inverse Gaussian as Bernstein distribution
#' parnamespcf("gigbd")

#' @export

parnamespcf = function(covtype,transform="none"){
  if (covtype=="expcov") parnames = "s"
  else if (covtype=="gammabd") parnames = c("tau","phi")
  else if (covtype=="invgammabd") parnames = c("tau","phi")
  else if (covtype=="gigbd") parnames = c("chi","psi","lambda")
  else stop("Unknown type of covariance function")
  if (transform=="lgcp") parnames = c("sigma",parnames)
  else if (transform=="icp") parnames = c("sigma","h",parnames)
  else if (transform=="pcpp") parnames = c("h",parnames)
  else if (transform!="none") stop("Unknown transform")
  return(parnames)
}
