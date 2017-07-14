## Graphical functions
# plotMaxMin function
#' @author Dario Azzimonti
#' @title Plot coordinate profiles
#' @name plotMaxMin
#' @description Plot coordinate profiles, for dimension up to 6.
#' @param allRes list containing the list \code{res} which contains the computed minima and maxima. The object returned by the function \code{getAllMaxMin}.
#' @param Design a d dimensional design corresponding to the points
#' @param threshold if not null plots the level
#' @param trueEvals if not null adds to each plot the data points and the observed value
#' @param km_model if not NULL plots the kriging model trend (TODO uncomment the lines).
#' @param ... additional parameters to be passed to the plot function
#' @return plots the sup and inf of the function for each dimension. If threshold is not NULL
# you can add km_model=NULL in signature and uncomment the lines() below to plot the trend.
# Still needs work on how to evaluate directly the formulas.
plotMaxMin<-function(allRes,Design=NULL,threshold=NULL,trueEvals=NULL,km_model=NULL,...){
  d<-ncol(allRes$res$min)
  if(is.null(Design)){
    if(!is.null(allRes$minima)){
      minR<-apply(allRes$minima,2,range)
      maxR<-apply(allRes$maxima,2,range)
      Design<-matrix(NA,ncol=d,nrow=100)
      for(i in seq(d)){
        Design[,i]<-seq(min(c(minR[,i],maxR[,i])),max(c(minR[,i],maxR[,i])),,100)
      }

    }else{
      warning("Design not selected, automatically chosen between 0 and 1")
      Design<-replicate(d,seq(0,1,,nrow(allRes$res$min)))
    }
  }

  mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
  oldpar<-par()
  # mar= c(5, 4, 4, 2) + 0.1
  par(mfrow=mfrows, mar = c(4, 3, 3, 1) + 0.1)
  for(coord in seq(d)){
    ylimTemp<-c(min(allRes$res$min[,coord]),max(allRes$res$max[,coord]))
    title_string<-paste("Coordinate",colnames(allRes$res$min)[coord])
    plot(Design[,coord],allRes$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
         xlab="x",ylab="f",...)
    lines(Design[,coord],allRes$res$max[,coord])
#    if(coord==3)
#      lines(Design[,coord],km_model@trend.coef[coord+1]*Design[,coord]^2,col=4)
#    lines(Design[,coord],km_model@trend.coef[coord+1]*Design[,coord],col=3)

    if(!is.null(threshold))
      abline(h = threshold,col=2)
    if(!is.null(trueEvals))
      points(trueEvals[,coord],trueEvals[,(d+1)])
  }
  par(mfrow=c(1,1))
}


# getChangePoints
#' @author Dario Azzimonti
#' @name getChangePoints
#' @title Coordinate profiles crossing points
#' @description Obtain the points where the coordinate profile extrema functions cross the threshold
#' @param allRes list containing the list \code{res} which contains the computed minima and maxima. The object returned by the function \code{getAllMaxMin}.
#' @param Design a d dimensional design corresponding to the points
#' @param threshold if not null plots the level
#' @return returns a list containing two lists with d components where
#' \itemize{
#' \item{alwaysEx:} each component is a numerical vector indicating the points \eqn{x_i} where \eqn{inf_{x^{-i}}f(x) >} \code{threshold};
#' \item{neverEx:} each component is a numerical vector indicating the points \eqn{x_i} where \eqn{sup_{x^{-i}}f(x) <} \code{threshold}.
#' }
getChangePoints <- function(threshold,Design=NULL,allRes){
  d<-ncol(allRes$res$min)
  chPointsAlwaysEx<-list()
  chPointsNeverEx<-list()
  if(is.null(Design))
    Design<-replicate(d,seq(0,1,,nrow(allRes$res$min)))

  alwaysExcursion <- allRes$res$min > threshold
  neverExcursion <- allRes$res$max < threshold

  for(i in seq(d)){
    if(!length(Design[alwaysExcursion[,i],i])){
      chPointsAlwaysEx[[i]]<-length(Design[alwaysExcursion[,i],i])
    }else{
      chPointsAlwaysEx[[i]]<-getSegments(Design[alwaysExcursion[,i],i])
    }

    if(!length(Design[neverExcursion[,i],i])){
      chPointsNeverEx[[i]]<-length(Design[neverExcursion[,i],i])
    }else{
      chPointsNeverEx[[i]]<-getSegments(Design[neverExcursion[,i],i])
    }

  }

  return(list(alwaysEx=chPointsAlwaysEx,neverEx=chPointsNeverEx))
}



# getSegments function
#' @author Dario Azzimonti
#' @name getSegments
#' @title getSegments
#' @description Auxiliary function for \code{getChangePoints}
#' @param y a vector
# @param Design a d dimensional design corresponding to the points
# @param threshold if not null plots the level
#' @return plots the sup and inf of the function for each dimension. If threshold is not NULL
getSegments = function(y){
  edges<-c(y[1])
  minDiffy<-min(unique(diff(y)))
  allJumps<-abs(diff(y)-minDiffy)>1e-15
  edges<-c(edges,y[allJumps])
  edges<-c(edges,y[c(FALSE,allJumps)],y[length(y)])
  return(sort(edges))
}
