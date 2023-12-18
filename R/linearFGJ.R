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
