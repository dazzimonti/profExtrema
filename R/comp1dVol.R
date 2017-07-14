#' @author Dario Azzimonti
#' @name compDm1Volume
#' @title d-1 dimensional volume of excursion
#' @description Computes the d-1 dimensional volume of excursion of \eqn{f} above \code{threshold}.
#' @param f the function of interest
#' @param threshold threshold of interest
#' @param d dimension of the input domain
#' @param opts a list containing the options for this function and the subfunction samplesFromCut
#' \itemize{
#' \item{\code{limits:}}{an optional list with the upper and lower limits of each dimension, if NULL then for each dimension limits are 0,1}
#' \item{\code{nSamples:}}{number of samples for MC integration in the d-1 dimensional space. If NULL \code{nSamples=200}.}
#' \item{\code{design:}}{a \eqn{d\times nDesign} design of experiments that defines the hyperplanes considered. If NULL \code{design[,coord]<-seq(0,1,,100)} for all coordinates.}
#' \item{\code{profileData:}}{Results on the profile maxima and minima. Useful to reduce the computations.}
#' }
#' @return a matrix with \code{d} columns and \code{length(opts$design)} rows. The element at (i,d) contains the volume of excursion on the d-1 dimensional hyperplane \eqn{x_d = design_i}.
compDm1Volume=function(f,threshold,d,opts=NULL){
  if(is.null(opts$limits)){
    opts$limits<-list(lower=rep(0,d),upper=rep(1,d))
  }
  if(is.null(opts$nSamples)){
    opts$nSamples<-200
  }
  flag=0
  if(!is.null(opts$profileData)){
    nDens<-length(opts$profileData$res$min[[1]])
    opts$design<-matrix(NA,ncol=d,nrow=nDens)
    flag=1
  }else{
    if(is.null(opts$design)){
      opts$design<-replicate(d,seq(0,1,,100))
    }
    nDens<-nrow(opts$design)
  }
  dens<-matrix(NA,ncol=d,nrow=nDens)

  for(coord in seq(d)){
    if(flag)
      opts$design[,coord]<-opts$profileData$minima[((coord-1)*nDens+1):(coord*nDens),coord]
    for(i in seq(nDens)){
      if(flag){
        curPoint<-opts$design[i,coord]
        if(opts$profileData$res$min[[coord]][i]>threshold){
          dens[i,coord]<-1
        }else if(opts$profileData$res$max[[coord]][i]<threshold){
          dens[i,coord]<-0
        }else{
          ss<-samplesFromCut(x=curPoint,
                             n = opts$nSamples,g = f,coord = coord,d = d,opts = opts)
          dens[i,coord]<-sum(ss$values>=threshold)/opts$nSamples
        }
      }else{
        curPoint<-opts$limits$lower[coord]+(opts$limits$upper[coord]-opts$limits$lower[coord])*opts$design[i,coord]
        ss<-samplesFromCut(x=curPoint,
                           n = opts$nSamples,g = f,coord = coord,d = d,opts = opts)
        dens[i,coord]<-sum(ss$values>=threshold)/opts$nSamples
      }

    }
  }
  if(is.null(opts$longReturn) || !opts$longReturn){
    return(dens)
  }else{
    return(list(volume=dens,design=opts$design))
  }
  return(dens)
}



#' @author Dario Azzimonti
#' @name samplesFromCut
#' @title Samples from 1d cut
#' @description Obtain samples from 1d cut
#' @param x to be documented
#' @param n to be documented
#' @param g to be documented
#' @param coord to be documented
#' @param d to be documented
#' @param opts a list containing the options
#' \itemize{
#' \item{\code{limits:}}{an optional list with the upper and lower limits of each dimension, if NULL then for each dimension limits are 0,1}
#' }
#' @return a list with \code{points} and \code{values}.
samplesFromCut<-function(x,n,g,coord,d,opts=NULL){
  if(is.null(opts)){
    opts$limits<-list(lower=rep(0,d),upper=rep(1,d))
  }
  pp<-sobol(n=n,dim=(d-1),scrambling = 1)
  #  pp<-matrix(runif(n*(d-1),min=opts$limits$lower[-coord],max=opts$limits$upper[-coord]),ncol=(d-1),nrow=n,byrow=T)
  if(coord==1){
    pp<-cbind(rep(x,n),pp)
  }
  if(coord==d){
    pp<-cbind(pp,rep(x,n))
  }
  if(coord>1 && coord<d){
    pp<-cbind(pp[,1:(coord-1)],rep(x,n),pp[,(coord):(d-1)])
  }
  rr<-apply(pp,1,g)
  return(list(points=pp,values=rr))
}
