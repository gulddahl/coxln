#' Find the relative boundary using resistance distance

#' @description
#' Find the relative boundary from a given point on a linear network
#' at a given resistance distance.

#' @details
#' This function finds the relative boundary from a point in a linear network,
#' and returns this as a point pattern on the same network. It uses the resistance
#' metric for calculating distances.

#' @param point A point pattern on a linear network consisting of a single point.
#' @param r The distance from the point to the relative boundary.
#' @param L For internal use.
#' @param ILM For internal use.
#' @param outpp If TRUE the a point pattern is returned; otherwise a dataframe is returned.
#' @return Either a point pattern of type lpp from spatstat or a dataframe.

#' @seealso [linearpcfR] for estimation of the pair correlation function
#' using the relative boundary for edge correction

#' @examples
#' # Relative boundary from point no. 30 in the spiders data
#' X = spatstat.data::spiders
#' point = X[30]
#' r = 180
#' rb = relativeboundaryR(point,r)
#' plot(point)
#' plot(rb,add=TRUE,cols=2)

#' # Relative boundary from point no. 10 in the chicago data
#' X = spatstat.data::chicago
#' point = X[10]
#' r = 100
#' rb = relativeboundaryR(point,r)
#' plot(point,use.marks=FALSE)
#' plot(rb,add=TRUE,cols=2)

#' @export

relativeboundaryR = function(point,r,L=as.linnet(point),ILM=solve(CalcLaplacianMatrix(L)),outpp=TRUE){
  segr = c(); tpr = c()
  j = as.numeric(point$data$seg)
  tpj = as.numeric(point$data$tp)
  l = lengths_psp(L$lines)
  Sigaiai = ILM[cbind(L$from,L$from)]
  Sigbibi = ILM[cbind(L$to,L$to)]
  Sigaibi = ILM[cbind(L$from,L$to)]
  Sigajai = as.vector(ILM[L$from[j],L$from])
  Sigajbi = as.vector(ILM[L$from[j],L$to])
  Sigaibj = as.vector(ILM[L$from,L$to[j]])
  Sigbibj = as.vector(ILM[L$to,L$to[j]])
  Sigajaj = rep(ILM[L$from[j],L$from[j]],nsegments(L))
  Sigbjbj = rep(ILM[L$to[j],L$to[j]],nsegments(L))
  Sigajbj = rep(ILM[L$from[j],L$to[j]],nsegments(L))
  Ai = (Sigaiai+Sigbibi-2*Sigaibi)/l^2-1/l
  Bijs = 1-2/l*(Sigaiai-Sigaibi-(1-tpj)*Sigajai+(1-tpj)*Sigajbi-tpj*Sigaibj+tpj*Sigbibj)
  Cijs = (1-tpj)^2*Sigajaj+tpj^2*Sigbjbj+2*tpj*(1-tpj)*Sigajbj+Sigaiai-2*(1-tpj)*Sigajai-2*tpj*Sigaibj+tpj*(1-tpj)*l[j]
  for (i in 1:nsegments(L)){
    if (i==j){ # case i=j
      if (Ai[i]>-1e-10){ # case Ai=0
        tpi1 = tpj+r/l[i]
        if (tpi1>=tpj&tpi1<=1) {tpr = c(tpr,tpi1); segr = c(segr,i)}
        tpi2 = tpj-r/l[i]
        if (tpi2>=0&tpi2<=tpj) {tpr = c(tpr,tpi2); segr = c(segr,i)}
      } else { # case Ai!=0
        d = 1+4*Ai[i]*r
        if (d<0) tpi1 = tpi2 = c()
        if (d==0) {tpi1 = -1/(2*Ai[i]*l[i])+tpj; tpi2 = 1/(2*Ai[i]*l[i])+tpj}
        if (d>0) {tpi1 = (-1+c(-1,1)*sqrt(d))/(2*Ai[i]*l[i])+tpj; tpi2 = (1+c(-1,1)*sqrt(d))/(2*Ai[i]*l[i])+tpj}
        tpr = c(tpr,tpi1[tpi1>=tpj&tpi1<=1]); segr = c(segr,rep(i,sum(tpi1>=tpj&tpi1<=1)))
        tpr = c(tpr,tpi2[tpi2>=0&tpi2<=tpj]); segr = c(segr,rep(i,sum(tpi2>=0&tpi2<=tpj)))
      }
    } else {  # case i!=j
      if (Ai[i]>-1e-10){  # case Ai=0
        tpi = (r-Cijs[i])/(Bijs[i]*l[i])
        tpr = c(tpr,tpi[tpi>=0&tpi<=1]); segr = c(segr,rep(i,sum(tpi>=0&tpi<=1)))
      } else {  # case Ai!=0
        d = Bijs[i]^2-4*Ai[i]*(Cijs[i]-r)
        if (d<0) tpi=c()
        if (d==0) tpi = -Bijs[i]/(2*Ai[i]*l[i])
        if (d>0) tpi = (-Bijs[i]+c(-1,1)*sqrt(d))/(2*Ai[i]*l[i])
        tpr = c(tpr,tpi[tpi>=0&tpi<=1]); segr = c(segr,rep(i,sum(tpi>=0&tpi<=1)))
      }
    }
  }
  if (outpp==TRUE) return(as.lpp(seg=segr,tp=tpr,L=L))
  else return(cbind(segr,tpr))
}


#' Calculate weights used in estimation of pcf

#' @description
#' Calculates the weights used in estimation of the pair correlation
#' function for a point pattern observed on a linear network using the
#' resistance metric.

#' @param u A point pattern on a linear network consisting of a single point.
#' @param r A distance.
#' @param L For internal use.
#' @param ILM For internal use.
#' @return A real number representing the weight.

#' @seealso [linearpcfR] for estimation of the pair correlation function
#' using the resistance metric

#' @export

wR = function(u,r,L=as.linnet(u),ILM=solve(CalcLaplacianMatrix(L))){
  rb = relativeboundaryR(u,r,L=L,ILM=ILM)
  J = abs(pairdistR(u,rb,derivative=TRUE,L=L,ILM=ILM))
  return(1/sum(1/J))
}


#' Estimate pcf using resistance metric

#' @description
#' Estimates the pair correlation function from a point pattern
#' on a linear network using the resistance metric.

#' @details
#' This function takes a point pattern on a linear network and estimates
#' the pair correlation function for a set of r-values by discretising the network.
#' The resistance metric is used for distances.

#' @param X A point pattern on a linear network.
#' @param r A vector of distances on which the estimate of the pair correlation
#' function is calculated.
#' @return The J function with class fv from spatstat.

#' @seealso [linearF], [linearG] and [linearJ] for other summary statistics

#' @examples
#' # Estimate pcf for spiders data
#' X = spatstat.data::spiders
#' r = seq(0,250,length.out=50)
#' pcf = linearpcfR(X,r)
#' plot(pcf)

#' @export

linearpcfR = function(X,r){
  L = as.linnet(X); n = npoints(X)
  ILM = solve(CalcLaplacianMatrix(L))
  pd = pairdistR(X)
  wRmat = matrix(0,n,n)
  for (i in 1:n){
    print(c(i,n))
    for (j in 1:n){
      if (i!=j) wRmat[i,j] = wR(u=X[i],r=pd[i,j],L=L,ILM=ILM)
    }
  }
  dens = density(as.vector(pd),bw="nrd0",weights = as.vector(wRmat),from=min(r),to=max(r),n=length(r),subdensity=TRUE,warnWbw=FALSE)
  pcfR = volume.linnet(L)/(n*(n-1)) * dens$y
  pcfR = fv(data.frame(r,pcfR),valu="pcfR",ylab="g(r)",yexp=quote(g(r)),fname="g")
  return(pcfR)
}

