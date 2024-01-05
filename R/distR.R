#' Extract diagonal from a matrix

#' @description
#' Returns the diagonal of a matrix, but unlike diag() it considers a
#' single number to be matrix.

#' @param x A square matrix.
#' @returns A vector with diagonal entries from x.

getdiag = function(x){
  diag(as.matrix(x))
}


#' Calculate the Laplacian matrix of a linear network

#' @description
#' Calculates the Laplacian matrix of a linear network from, and adds a +1
#' to the first entry in the matrix to make it invertible.

#' @param L A linear network of class linnet from the spatstat package.
#' @returns A matrix.

#' @seealso [distR] and [pairdistR] for calculation of the resistance
#' metric using the Laplacian matrix

#' @examples
#' CalcLaplacianMatrix(spatstat.data::simplenet)
#'
#' L = as.linnet(spatstat.data::spiders)
#' CalcLaplacianMatrix(L)
#'
#' L = as.linnet(spatstat.data::dendrite,sparse=FALSE)
#' CalcLaplacianMatrix(L)

#' @export

CalcLaplacianMatrix = function(L){
  if (L$sparse==TRUE) stop("Error: L has sparse representation (may be fixed with L = as.linnet(L,sparse=FALSE)).")
  Delta = -L$m*(1/L$dpath)
  diag(Delta) = 0
  diag(Delta) = -apply(Delta,1,sum)+c(1,rep(0,dim(Delta)[1]-1))
  return(Delta)
}


#' Resistance distance between two points in a linear network

#' @description
#' Calculates the resistance distance between two points u and v in a
#' linear network.

#' @param segu The line segment number of point u.
#' @param tu The position of u on the line segment (a value in [0,1])
#' @param segv The line segment number of point v.
#' @param tv The position of v on the line segment (a value in [0,1])
#' @param L A linear network of class linnet from the spatstat package.
#' @param ILM For internal use.
#' @returns The resistance distance between u and v.

#' @seealso [pairdistR] for calculation of the resistance
#' metric on one or two point pattern on a linear network,
#' [CalcLaplacianMatrix] for calculation of the Laplacian matrix

#' @examples
#' # Distance between midpoints in segments number 1 and 2
#' distR(1,0.5,2,0.5,spatstat.data::simplenet)

#' @export

distR = function(segu,tu,segv,tv,L,ILM = solve(CalcLaplacianMatrix(L))){
  Sigajaj = ILM[L$from[segu],L$from[segu]]
  Sigajbj = ILM[L$from[segu],L$to[segu]]
  Sigbjbj = ILM[L$to[segu],L$to[segu]]
  Sigaiai = ILM[L$from[segv],L$from[segv]]
  Sigaibi = ILM[L$from[segv],L$to[segv]]
  Sigbibi = ILM[L$to[segv],L$to[segv]]
  Sigajai = ILM[L$from[segu],L$from[segv]]
  Sigajbi = ILM[L$from[segu],L$to[segv]]
  Sigbjai = ILM[L$to[segu],L$from[segv]]
  Sigbjbi = ILM[L$to[segu],L$to[segv]]
  li = lengths_psp(L$lines)[segv]
  lj = lengths_psp(L$lines)[segu]
  dist = ((1-tu)^2 * Sigajaj + tu^2 * Sigbjbj + 2*tu*(1-tu) * Sigajbj
          + (1-tv)^2 * Sigaiai + tv^2 * Sigbibi + 2*tv*(1-tv) * Sigaibi
          - 2*(1-tu)*(1-tv) * Sigajai - 2*(1-tu)*tv * Sigajbi
          - 2*tu*(1-tv) * Sigbjai - 2*tu*tv * Sigbjbi
          + tv*(1-tv)*li + tu*(1-tu)*lj
          - 2*(segu==segv)*min(tu*(1-tv),tv*(1-tu))*li)
  return(dist)
}


#' Pairwise distances of point pattern on linear network using the resistance metric

#' @description
#' Calculates all pairwise distances using the resistance distance between all points in
#' a point pattern or two patterns on a linear network. Can also calculate the derivative
#' (i.e. Jacobian) instead.

#' @param X If a single point pattern is provided all pairwise distances (using the
#' resistance distance) between points in X are found.
#' @param X2 If a second point pattern is provided, all pairwise distances will be found
#' between points in X and X2 are found.
#' @param derivative If TRUE then the derivative of the resistance distance d(u,v) wrt. v
#' is found instead of the resistance distance itself. Default is FALSE.
#' @param L Internal use.
#' @param ILM Internal use.
#' @returns A matrix with all the distances or derivatives of the distances.

#' @seealso [distR] for calculation of the resistance
#' metric between two points in a linear network,
#' [CalcLaplacianMatrix] for calculation of the Laplacian matrix

#' @examples
#' # All pairwise distances in the spiders dataset
#' X = spatstat.data::spiders
#' pairdistR(X)
#'
#' # Sparse representations in the linear network does not work
#' X = spatstat.data::dendrite
#' X = as.lpp(X,sparse=FALSE)
#' pairdistR(X)
#'
#' # All pairwise distances between two Poisson processes on simplenet
#' L = spatstat.data::simplenet
#' X = rpoislpp(5,L)
#' Y = rpoislpp(5,L)
#' pairdistR(X,Y)
#'
#' # The derivative of the resistance distance
#' X = spatstat.data::spiders
#' pairdistR(X,derivative=TRUE)

#' @export

pairdistR = function(X,X2=NULL,derivative=FALSE,L=as.linnet(X),ILM=solve(CalcLaplacianMatrix(L))){
  if (is.null(X2)){
    segu = segv = as.integer(X$data$seg); tu = tv = as.numeric(X$data$tp)
  } else {
    segu = as.integer(X$data$seg); tu = as.numeric(X$data$tp)
    segv = as.integer(X2$data$seg); tv = as.numeric(X2$data$tp)
  }
  onesu = rep(1,length(segu))
  onesv = rep(1,length(segv))
  Sigajaj = outer(getdiag(ILM)[L$from[segu]],onesv)
  Sigajbj = outer(getdiag(ILM[L$from[segu],L$to[segu]]),onesv)
  Sigbjbj = outer(getdiag(ILM)[L$to[segu]],onesv)
  Sigaiai = outer(onesu,getdiag(ILM)[L$from[segv]])
  Sigaibi = outer(onesu,getdiag(ILM[L$from[segv],L$to[segv]]))
  Sigbibi = outer(onesu,getdiag(ILM)[L$to[segv]])
  Sigajai = ILM[L$from[segu],L$from[segv]]
  Sigajbi = ILM[L$from[segu],L$to[segv]]
  Sigbjai = ILM[L$to[segu],L$from[segv]]
  Sigbjbi = ILM[L$to[segu],L$to[segv]]
  li = outer(onesu,lengths_psp(L$lines)[segv])
  liv = lengths_psp(L$lines)[segv]
  lj = outer(lengths_psp(L$lines)[segu],onesv)
  tum = outer(tu,onesv)
  tvm = outer(onesu,tv)
  if (derivative==FALSE){
    dist = ((1-tum)^2 * Sigajaj + tum^2 * Sigbjbj + 2*tum*(1-tum) * Sigajbj
            + (1-tvm)^2 * Sigaiai + tvm^2 * Sigbibi + 2*tvm*(1-tvm) * Sigaibi
            - 2*(1-tum)*(1-tvm) * Sigajai - 2*(1-tum)*tvm * Sigajbi
            - 2*tum*(1-tvm) * Sigbjai - 2*tum*tvm * Sigbjbi
            + tvm*(1-tvm)*li + tum*(1-tum)*lj
            - 2*outer(segu,segv,"==")*pmin(outer(tu,1-tv),outer(1-tu,tv))*outer(onesu,liv))
    return(dist*(dist>0))   ### removing negative values due to numerical rounding
  } else {
    J = (- 2*(1-tvm)*Sigaiai/li + 2*tvm*Sigbibi/li + 2*(1-2*tvm)*Sigaibi/li
         + 2*(1-tum)*Sigajai/li - 2*(1-tum)*Sigajbi/li + 2*tum*Sigbjai/li - 2*tum*Sigbjbi/li
         + (1-2*tvm) + 2*tum*outer(segu,segv,"==")*(outer(tu,1-tv)<outer(1-tu,tv))
         - 2*(1-tum)*outer(segu,segv,"==")*(outer(tu,1-tv)>outer(1-tu,tv)))
    return(J)
  }
}

