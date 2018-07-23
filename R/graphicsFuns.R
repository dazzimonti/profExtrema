## Graphical functions
# plotMaxMin function
#' @author Dario Azzimonti
#' @title Plot coordinate profiles
#' @name plotMaxMin
#' @description Plot coordinate profiles, for dimension up to 6.
#' @param allRes list containing the list \code{res} which contains the computed minima and maxima. The object returned by the function \code{getAllMaxMin}.
#' @param Design a d dimensional design corresponding to the points
#' @param threshold if not \code{NULL} plots the level
#' @param changes boolean, if not \code{FALSE} plots the points where profile extrema take values near the threshold.
#' @param trueEvals if not \code{NULL} adds to each plot the data points and the observed value
#' @param ... additional parameters to be passed to the plot function
#' @return plots the sup and inf of the function for each dimension. If threshold is not NULL
#' @export
# you can add km_model=NULL in signature and uncomment the lines() below to plot the trend.
# Still needs work on how to evaluate directly the formulas.
plotMaxMin<-function(allRes,Design=NULL,threshold=NULL,changes=FALSE,trueEvals=NULL,...){
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

  if(changes){
    pp_change<-getChangePoints(threshold = threshold,Design = Design,allRes=allRes)
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

    if(!is.null(threshold))
      abline(h = threshold,col=2)
    if(changes){
      for(cc in seq(names(pp_change$alwaysEx))){
        abline(v = pp_change$alwaysEx[[cc]][[coord]],col=4)
        abline(v = pp_change$neverEx[[cc]][[coord]],col="darkgreen")
      }
    }
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
#' @export
getChangePoints <- function(threshold,Design=NULL,allRes){
  d<-ncol(allRes$res$min)
  if(is.null(Design))
    Design<-replicate(d,seq(0,1,,nrow(allRes$res$min)))

  # adapt the function to many thresholds
  num_t<-length(threshold)
  alwaysExcursion<-array(NA,dim = c(nrow(Design),d,num_t))
  neverExcursion<-array(NA,dim = c(nrow(Design),d,num_t))

  chPointsAlwaysEx<-replicate(n = num_t,list())
  names(chPointsAlwaysEx)<-round(threshold)
  chPointsNeverEx<-replicate(n = num_t,list())
  names(chPointsNeverEx)<-round(threshold)

  for(j in seq(num_t)){
    alwaysExcursion[,,j] <- allRes$res$min > threshold[j]
    neverExcursion[,,j] <- allRes$res$max < threshold[j]

    for(i in seq(d)){
      if(!length(Design[alwaysExcursion[,i,j],i])){
        chPointsAlwaysEx[[j]][[i]]<-min(Design[,i])
      }else{
        chPointsAlwaysEx[[j]][[i]]<-getSegments(Design[alwaysExcursion[,i,j],i])
      }

      if(!length(Design[neverExcursion[,i,j],i])){
        chPointsNeverEx[[j]][[i]]<-min(Design[,i])
      }else{
        chPointsNeverEx[[j]][[i]]<-getSegments(Design[neverExcursion[,i,j],i])
      }

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
#' @export
getSegments = function(y){
  edges<-c(y[1])
  minDiffy<-min(unique(diff(y)))
  allJumps<-abs(diff(y)-minDiffy)>1e-15
  edges<-c(edges,y[allJumps])
  edges<-c(edges,y[c(FALSE,allJumps)],y[length(y)])
  return(sort(edges))
}



# plotOblique function
#' @author Dario Azzimonti
#' @name plotOblique
#' @title plotOblique
#' @description Auxiliary function for 2d plotting of excluded regions
#' @param changePoints Numerical vector with the change points (usually if \code{cp=getChangePoints(...)}, then this is cc$alwaysEx[[1]][[1]] for example)
#' @param direction The Psi vector used for the direction
#' @param ... parameters to be passed to abline
#' @return adds to the current plot the lines \eqn{x} s.t. \code{direction}^T \eqn{x} = \code{changePoints[i]} for all i
#' @export
plotOblique<-function(changePoints,direction,...){

  av<- changePoints/direction[2]
  bv<- -direction[1]/direction[2]*(changePoints!=0)

  for(i in seq(length(changePoints))){
    if(direction[2]==0){
      abline(v=changePoints[i]/direction[1],...)
    }else{
      if(av[i]==0){
        next
      }else{
        abline(a=av[i],b=bv[i],...)
      }
    }

  }
}


# plotOneBivProfile function
#' @author Dario Azzimonti
#' @title Plot bivariate profiles
#' @name plotOneBivProfile
#' @description Plots the bivariate profiles stored in \code{allRes} for each Psi in \code{allPsi}.
#' @param allRes list containing the list \code{res} which contains the computed minima and maxima. The object returned by the function \code{getProfileExtrema}.
#' @param allPsi a list containing the matrices Psi (dim \eqn{2xd}) for which to compute the profile extrema
#' @param Design a matrix of dimension \eqn{(2d)x numPsi} encoding the first (\code{Design[1:d,]}) and the second ((\code{Design[(d+1):(2*d),]})) axis values.
#' @param threshold if not \code{NULL} plots the level as a contour.
#' @param trueEvals if not \code{NULL} adds to each plot the data points and the observed value
#' @param main_addendum additional string to add to image title. Default is empty string.
#' @param ... additional parameters to be passed to the plot function
#' @return plots the 2d maps of the profile sup and inf in \code{allRes} for each Psi in \code{allPsi}. If threshold is not NULL also contours the threshold level.
#' @seealso plotBivariateProfiles
#' @export
plotOneBivProfile<-function(allRes,allPsi,Design=NULL,threshold=NULL,trueEvals=NULL,main_addendum="",...){


  num_Psi<-length(allPsi)



  if(is.null(Design)){
      warning("Design not selected, auto choice not implemented")
      return(NULL)
  }

  dd_eta<-length(Design[,1])/2

  bkg_col <- gray.colors(20)

  for(i in seq(num_Psi)){

    zzBackground_max <- matrix(allRes$res$max[[i]],ncol=dd_eta)
    zzBackground_min <- matrix(allRes$res$min[[i]],ncol=dd_eta)

#    par(mfrow=c(1,2))
    layout(matrix(1:4,nrow=1),widths=c(0.45,0.05,0.45,0.05))
    par(mar=c(4,4,4,0.2) + 0.2)
    image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackground_max,col=bkg_col,
          main=bquote(P[Psi[.(i)]]^sup ~ "f " ~ .(main_addendum)),...)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackground_max,nlevels = 10,add=T)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackground_max,levels = threshold,add=T,col=2,lwd=1.4)

    legend_info<-hist(zzBackground_max,plot=FALSE,breaks=20)
    which_col_leg <- which(legend_info$counts>0)
    par(mar=c(5,0.2,5,2.8))
    image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
    axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))

    par(mar=c(4,4,4,0.2) + 0.2)
    image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackground_min,col=bkg_col,main=bquote(P[Psi[.(i)]]^inf ~ "f" ~ .(main_addendum)),...)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackground_min,nlevels = 10,add=T)
    contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=zzBackground_min,levels = threshold,add=T,col=2,lwd=1.4)

    legend_info<-hist(zzBackground_min,plot=FALSE,breaks=20)
    which_col_leg <- which(legend_info$counts>0)
    par(mar=c(5,0.2,5,2.8))
    image(y=1:length(which_col_leg),z=t(1:length(which_col_leg)), col=bkg_col[which_col_leg], axes=FALSE)#, main="Slope", cex.main=.8)
    axis(4,cex.axis=0.8,labels = legend_info$breaks[which_col_leg],at=1:length(which_col_leg),...)#,mgp=c(0,.5,0))
  }
#  par(mfrow=c(1,1))


}


