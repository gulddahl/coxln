#' Find end points in network

#' @description
#' Finds all the end points in a network, i.e. all the vertices
#' of degree one in the graph underlying the network.

#' @details
#' Given a linear network, this function returns a point pattern
#' on the same network giving all the end points in the network.

#' @param L A linear network.
#' @returns A point pattern on a linear network with class lpp
#' from spatstat.

#' @examples
#' # Find the end points from simplenet from spatstat
#' ep = findendpoints(spatstat.data::simplenet)
#' plot(ep)

#' @export

findendpoints = function(L){
  x = sort(c(L$from,L$to))
  x = unique(x*duplicated(x))
  x = x[x>0]
  epi = (1:npoints(L$vertices))[-x]
  ep = as.lpp(L$vertices[epi],L=L)
  return(ep)
}


#' Removes points within distance r of end points

#' @description
#' Removes all points within distance r from the end points of
#' a point pattern defined on a linear network.

#' @details
#' This function takes a point pattern X on a linear network, and
#' removes all points within distance r from the end points of the
#' network, i.e. all the vertices of degree one in the graph
#' underlying the network. The distances can be measured either by
#' geodesic distance or resistance distance.

#' @param X A point pattern on a linear network.
#' @param r The distance to the end points, within which points should be
#' removed.
#' @param metric The metric used to measure distances. Should be either "G"
#' for geodesic distance, or "R" for resistance distance.
#' @param saveind For internal use.
#' @return A point pattern on a linear network with class lpp from spatstat.

#' @examples
#' # Removes points from the dendrite dataset using geodesic metric
#' X = spatstat.data::dendrite
#' r = 20
#' Xm = minussampling(X,r,metric="G")
#' plot(X)
#' plot(Xm, add=TRUE, cols=3)
#'
#' # Removes points from the spiders dataset using resistance metric
#' X = spatstat.data::spiders
#' r = 200
#' Xm = minussampling(X,r,metric="R")
#' plot(X)
#' plot(Xm, add=TRUE, cols=3)

#' @export

minussampling = function(X,r,metric="G",saveind=FALSE){
  ep = findendpoints(as.linnet(X))
  if (metric=="G") {ind = distfun.lpp(ep)(X)>r}  # geodesic metric
  else {ind = apply(pairdistR(X,ep),1,min)>r}   # resistance metric
  Xm = X[ind]
  if (saveind==TRUE) attr(Xm, "ind") = ind
  return(Xm)
}
