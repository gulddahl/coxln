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


#' Estimate F function on linear network

#' @description
#' Estimates the F function for a point pattern observed on a linear network.

#' @details
#' This function takes a point pattern on a linear network and estimates
#' the F function for a set of r-values by discretising the network. Both the
#' geodesic and the resistance metric can be used for distances

#' @param X A point pattern on a linear network.
#' @param r A vector of distances on which the estimate of F(r) is calculated.
#' @param ppul The number of points per unit length of the linear network used
#' for approximating F. This is rounded down to an integer for each line segment.
#' Note that the vertices in the graph underlying the linear network are
#' automatically added, and are not counted in ppul.
#' @param metric The metric used to measure distances. Should be either "G"
#' for geodesic distance, or "R" for resistance distance.
#' @return The F function with class fv from spatstat.

#' @examples
#' # Estimate F function for spiders data
#' X = spatstat.data::spiders
#' r = seq(0,250,length.out=50)
#' ppul = 0.02
#' F = linearF(X,r,ppul,metric="G")
#' plot(F)
#'
#' # Estimate F function for dendrite data
#' X = spatstat.data::dendrite
#' r = seq(0,5,length.out=50)
#' ppul = 0.5
#' F = linearF(X,r,ppul,metric="G")
#' plot(F)

#' @export

linearF = function(X,r,ppul,metric="G"){
  nr = length(r)
  Fest = rep(0,nr)
  pos = makepos(as.linnet(X),ppul)
  for (i in 1:nr){
    posm = minussampling(pos,r[i],metric=metric)
    if (metric=="G") Fest[i] = 1-npoints(posm[distfun.lpp(X)(posm)>r[i]])/npoints(posm)
    else Fest[i] = 1-sum(apply(pairdistR(posm,X),1,min)>r[i])/npoints(posm)
  }
  Fest = fv(data.frame(r,Fest),valu="Fest",ylab="F(r)",yexp=quote(F(r)),fname="F")
  return(Fest)
}


#' Estimate G function on linear network

#' @description
#' Estimates the G function for a point pattern observed on a linear network.

#' @details
#' This function takes a point pattern on a linear network and estimates
#' the G function for a set of r-values. Both the geodesic and the
#' resistance metric can be used for distances

#' @param X A point pattern on a linear network.
#' @param r A vector of distances on which the estimate of G(r) is calculated.
#' @param metric The metric used to measure distances. Should be either "G"
#' for geodesic distance, or "R" for resistance distance.
#' @return The G function with class fv from spatstat.

#' @examples
#' # Estimate G function for spiders data
#' X = spatstat.data::spiders
#' r = seq(0,250,length.out=50)
#' G = linearG(X,r,metric="G")
#' plot(G)
#'
#' # Estimate G function for dendrite data
#' X = spatstat.data::dendrite
#' r = seq(0,5,length.out=50)
#' G = linearG(X,r,ppul,metric="G")
#' plot(G)

#' @export

linearG = function(X,r,metric="G"){
  nr = length(r)
  Gest = rep(1,nr)
  Xm = X
  if (metric=="G") pd = pairdist.lpp(Xm) else pd = pairdistR(Xm)
  for (i in 1:nr){
    Xm = minussampling(Xm,r[i],metric=metric,saveind=TRUE)
    ind = attr(Xm,"ind")
    pd = pd[ind,ind]
    if (npoints(Xm)<2) {
      Gest[i] = 1
      break  # stop for loop if there are less than two points left
    }
    else{
      diag(pd) = Inf
      Gest[i] = 1-sum(apply(pd,1,min)>r[i])/npoints(Xm)
    }
  }
  Gest = fv(data.frame(r,Gest),valu="Gest",ylab="G(r)",yexp=quote(G(r)),fname="G")
  return(Gest)
}

#' Estimate J function on linear network

#' @description
#' Estimates the J function for a point pattern observed on a linear network.

#' @details
#' This function takes a point pattern on a linear network and estimates
#' the J function for a set of r-values by discretising the network. Both the
#' geodesic and the resistance metric can be used for distances

#' @param X A point pattern on a linear network.
#' @param r A vector of distances on which the estimate of J(r) is calculated.
#' @param ppul The number of points per unit length of the linear network used
#' for approximating the F function used in the J-function.
#' This is rounded down to an integer for each line segment.
#' Note that the vertices in the graph underlying the linear network are
#' automatically added, and are not counted in ppul.
#' @param metric The metric used to measure distances. Should be either "G"
#' for geodesic distance, or "R" for resistance distance.
#' @return The J function with class fv from spatstat.

#' @examples
#' # Estimate J function for spiders data
#' X = spatstat.data::spiders
#' r = seq(0,250,length.out=50)
#' ppul = 0.02
#' J = linearJ(X,r,ppul,metric="G")
#' plot(J)
#'
#' # Estimate J function for dendrite data
#' X = spatstat.data::dendrite
#' r = seq(0,5,length.out=50)
#' ppul = 0.5
#' J = linearJ(X,r,ppul,metric="G")
#' plot(J)

#' @export

linearJ = function(X,r,ppul,metric="G"){
  F = linearF(X,r,ppul,metric=metric)
  G = linearG(X,r,metric=metric)
  Jest = (1-G$Gest)/(1-F$Fest)
  Jest = fv(data.frame(r,Jest),valu="Jest",ylab="J(r)",yexp=quote(J(r)),fname="J")
  return(Jest)
}

