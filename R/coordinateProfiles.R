# coordinateProfiles function
#' @author Dario Azzimonti
#' @title Coordinate profiles starting from a kriging model
#' @description The function coordinateProfiles computes the profile extrema functions for the posterior mean of a Gaussian process and its confidence bounds
#' @param kmModel a \link[DiceKriging]{km} model
#' @param threshold the threshold of interest
#' @param options_full an optional list of options for getAllMaxMin, see \link{getAllMaxMin} for details.
#' @param options_approx an optional list of options for approxMaxMin, see \link{approxMaxMin} for details.
#' @param uq_computations boolean, if TRUE the uq computations for the profile mean are computed.
#' @param plot_level an integer to select the plots to return (0=no plots, 1=basic plots, 2= all plots)
#' @param plot_options an optional list of parameters for plots. Currently available options
#' \itemize{
#' \item{\code{save:}}{boolean, if TRUE saves the plots in \code{folderPlots}}
#' \item{\code{folderPlots:}}{a string containing the destination folder for plots, if \code{save==TRUE} default is \code{./}}
#' \item{\code{titleProf:}}{a string containing the title for the coordinate profile plots}
#' \item{\code{title2d:}}{a string containing the title for the 2d plots (if the input is 2d)}
#' \item{\code{design:}}{a \eqn{dxr} matrix where \eqn{d} is the input dimension and \eqn{r} is the size of the discretization for plots at each dimension}
#' \item{\code{id_save:}}{a string to be added to the plot file names, useful for serial computations on HPC.}
#' \item{\code{qq_fill:}}{if TRUE it fills the region between the first 2 quantiles in \code{quantiles_uq}.}
#' }
#' @param return_level an integer to select the amount of details returned
#' @param ... additional parameters to be passed to \link{coordProf_UQ}.
#' @return If return_level=1 a list containing
#' \itemize{
#' \item{\code{profMean_full:}}{the results of \code{getAllMaxMin} for the posterior mean}
#' \item{\code{profMean_approx:}}{the results of \code{approxMaxMin} for the posterior mean}
#' \item{\code{res_UQ:}}{the results of \code{coordProf_UQ} for the posterior mean}
#' }
#' if return_level=2 the same list as above but also including
#' \itemize{
#' \item{\code{abs_err:}}{the vector of maximum absolute approximation errors for the profile inf /sup on posterior mean for the chosen approximation}
#' \item{\code{times:} }{ a list containing
#' 	\itemize{
#' 	\item{\code{full:}}{computational time for the full computation of profile extrema}
#' 	\item{\code{approx:}}{computational time for the approximate computation of profile extrema}
#' 	}}
#' }
#' @export
coordinateProfiles = function(kmModel,threshold,options_full=NULL,options_approx=NULL,uq_computations=FALSE,plot_level=0,plot_options=NULL,return_level=1,...){

  # Options setup
  if(is.null(options_full)){
    options_full<-list(multistart=8,heavyReturn=TRUE)
  }

  if(is.null(options_approx)){
    init_des<-maximinLHS(10,2)
    options_approx<- list(multistart=8,heavyReturn=TRUE,initDesign=init_des,fullDesignSize=100,smoother="1order")
  }

  if(is.null(plot_options)){
    if(plot_level>0){
      plot_options<-list(save=F)
    }

  }else{
    if(plot_options$save==TRUE && is.null(plot_options$folderPlots))
      plot_options$folderPlots <- './'

    if(is.null(plot_options$titleProf))
      plot_options$titleProf<-"Coordinate profiles"

    if(is.null(plot_options$title2d))
      plot_options$title2d<-"Posterior mean"

    if(is.null(plot_options$design)){
      plot_options$design<-matrix(NA,ncol=kmModel@d,nrow=100)
      for(i in seq(kmModel@d)){
        plot_options$design[,i]<-seq(0,1,,100)
      }
    }
    # Useful for serial plots in HPC
    if(is.null(plot_options$id_save)){
      plot_options$id_save<-""
    }

  }

  # save the dimension of the input
  d<-kmModel@d

  # Let us define the posterior mean function and gradient in a 'optimizer-friendly' way
  g <-function(y){
    y<-matrix(y,ncol=d)
    return(predict.km(object = kmModel,newdata = y,type="UK",light.return = TRUE,checkNames = FALSE)$mean)
  }
  gprime <- function(y){
    y<-matrix(y,ncol=d)
    return(gradKm_dnewdata(object = kmModel,newdata = y,type = "UK",light.return = TRUE)$mean)
  }

  # Compute profile extrema with full optim
  timeIn<-get_nanotime()
  allRes_full<-getAllMaxMin(f = g,fprime = gprime,d = d,options = options_full)
  timeMean_full<-(get_nanotime()-timeIn)*1e-9


  # Obtain the threshold points selected by profile extrema
  changePP<-getChangePoints(threshold = threshold,allRes = allRes_full)


  # Plot posterior mean, ONLY for dimension==2
  if(plot_level>=2 && d==2){
    # since dimension==2 we can plot the posterior mean
    newdata<-expand.grid(seq(0,1,,100),seq(0,1,,100))
    colnames(newdata)<-colnames(kmModel@X)
    pred2d<-predict.km(kmModel,newdata = newdata,type = "UK",light.return = TRUE,se.compute = FALSE)

    # Plot the posterior mean and visualize the actual excursion set and the regions of no-excursion according to the profile extrema functions.
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"2dpostMean_1",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    par(mar = c(5, 5, 4, 2) + 0.1)
    image(matrix(pred2d$mean,nrow = 100),col=gray.colors(20), main=plot_options$title2d,xlab = colnames(kmModel@X)[1],ylab= colnames(kmModel@X)[2],
          cex.main=3,cex.axis=1.8,cex.lab=2.8)
    contour(matrix(pred2d$mean,nrow = 100),add=T,nlevels = 10,lwd=1.5,labcex=1.2)
    contour(matrix(pred2d$mean,nrow = 100),add=T,levels = threshold,col=2,lwd=3,labcex=1.5)
    abline(v = changePP$neverEx[[1]],col=3,lwd=2.5)
    abline(h = changePP$neverEx[[2]],col=3,lwd=2.5)
    if(plot_options$save)
      dev.off()
  }

  # Plot the profile extrema functions
  if(plot_level>0){
    colnames(allRes_full$res$min)<-colnames(kmModel@X)
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"profMean_full",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    plotMaxMin(allRes = allRes_full,threshold = threshold)
    if(plot_options$save)
      dev.off()
  }

  # Compute profile extrema with approximation
  timeIn<-get_nanotime()
  allRes_approx<-approxMaxMin(f = g,fprime = gprime,d = kmModel@d,opts = options_approx)
  timeMean_approx<-(get_nanotime()-timeIn)*1e-9


  if(plot_level>0){
    oldpar<-par()
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"profMean_comparison",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
    par(mfrow=mfrows, mar = c(4, 5, 3, 1) + 0.1)
    for(coord in seq(d)){
      ylimTemp<-c(min(allRes_full$res$min[,coord]),max(allRes_full$res$max[,coord]))
      title_string<-paste("Coordinate",colnames(allRes_full$res$min)[coord])
      plot(plot_options$design[,coord],allRes_full$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
           xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      lines(plot_options$design[,coord],allRes_full$res$max[,coord],lty=1)
      lines(plot_options$design[,coord],allRes_approx$res$max[,coord],lty=3,col=4,lwd=2)
      points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$max[,coord],col=4)
      lines(plot_options$design[,coord],allRes_approx$res$min[,coord],lty=4,col=4,lwd=2)
      points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$min[,coord],col=4)
      abline(v=changePP$neverEx[[coord]],col=3,lwd=2.5)
      abline(h = threshold,col=2,lwd=2)
      legend("bottomleft",c(as.expression(substitute(paste(P[coord]^sup,f," (full)"),list(coord=coord))),
                            as.expression(substitute(paste(P[coord]^inf,f," (full)"),list(coord=coord))),
                            as.expression(substitute(paste(P[coord]^sup,f," (1order)"),list(coord=coord))),
                            as.expression(substitute(paste(P[coord]^inf,f," (1order)"),list(coord=coord))),
                            "threshold"),
             lty=c(1,1,3,4,1),col=c(1,1,4,4,2),cex=1.8,lwd=c(2,2,2,2,2))
    }
    if(plot_options$save)
      dev.off()
  #  par(oldpar)
  }

  # save absolute errors with approximation
  abs_err_min<-max(abs(allRes_approx$res$min-allRes_full$res$min)/allRes_full$res$min)
  abs_err_max<-max(abs(allRes_approx$res$max-allRes_full$res$max)/allRes_full$res$max)
  abs_err<-c(abs_err_min,abs_err_max)

  ## UQ part #
  res_UQ<-NULL
  if(uq_computations){
    res_UQ<-coordProf_UQ(kmModel=kmModel,threshold=threshold,allResMean=allRes_full,
                         options_approx=options_approx,plot_level=max(0,plot_level-1),
                         plot_options=plot_options,return_level=max(1,return_level-1),...)

  }
  results<-list(profMean_full=allRes_full,profMean_approx=allRes_approx,res_UQ=res_UQ)
  if(return_level==1){
    return(results)
  }else{
    times<-list(full=timeMean_full,approx=timeMean_approx)
    results$more<-list(abs_err=abs_err,times=times)
    return(results)
  }

#  save(allRes_approx,allRes_full,Design,profInfs,profSups,kmModel,m_dist,some.simu,file=paste(plot_options$folderPlots,"Res_2dBRGM.RData",sep=""))
}
