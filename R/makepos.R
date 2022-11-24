#' Make a regular grid of points on a linear network

#' @description
#' makepos() makes a grid of points on a linear network
#' by placing points on the vertices and points regularly
#' spaced on each line segment.

#' @param L A linear network of class linnet from the spatstat package.
#' @param ppul Number of points per unit length on each line segment.
#' @param duplicate If duplicate is TRUE each vertex will contain
#'   duplicate points (one for each line segment extending
#'   from the vertex); if duplicate is FALSE, there is only one
#'   point on each vertex.
#' @returns A point pattern of class lpp from the spatstat package.

#' @examples
#' makepos(spatstat.data::simplenet,10)

#' @export

makepos = function(L,ppul,duplicate=FALSE){
  tp = seg = c()
  l = lengths_psp(L$lines)
  n = floor(ppul*l)
  for (i in 1:nsegments(L)){
    tpi = 0:(n[i]+1)/(n[i]+1)
    tp = c(tp,tpi)
    seg = c(seg,rep(i,n[i]+2))
  }
  pos = as.lpp(tp=tp,seg=seg,L=L)
  if (duplicate==TRUE) {
    return(pos)
  } else {
    dist = pairdist.lpp(pos)
    dist[upper.tri(dist,diag=TRUE)]=Inf
    pos2 = pos[apply(dist,1,min)>0]
    return(pos2)
  }
}

