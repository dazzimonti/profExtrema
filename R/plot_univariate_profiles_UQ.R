# Function to plot the univariate profile extrema functions with UQ
#' @name plot_univariate_profiles_UQ
#' @description Function to plot the univariate profile extrema functions with UQ
#' @title Univariate profile extrema with UQ
#' @param objectUQ an object returned by \link{coordProf_UQ} or the object saved in \code{obj$res_UQ}, if \code{obj} is the object returned by \link{coordinateProfiles}.
#' @param plot_options a list containing the same elements as the one passed to \link{coordinateProfiles}
#' @param nsims number of simulations
#' @param threshold threshold of interest
#' @param nameFile the central name of the plot file
#' @param quantiles_uq a vector containing the quantiles to be computed
#' @param profMean the profile coordinate extrema functions for the mean. It is saved in \code{obj$profMean_full} or \code{obj$profMean_approx} if \code{obj} is an object returned by \link{coordinateProfiles}.
#' @param typeProf a string to choose with type of profile extrema for simulations to plot \itemize{
#' \item{"approx"}{plots only the approximate profile extrema for simulations}
#' \item{"full"}{plots only the full profile extrema for simulations}
#' \item{"both"}{plots both the approximate and full profile extrema for simulations}
#' }
#' @return Plots either to the default graphical device or to pdf (according to the options passed in \code{plot_options})
#' @export
plot_univariate_profiles_UQ<-function(objectUQ,plot_options,nsims,threshold,nameFile="prof_UQ",quantiles_uq=c(0.05,0.95),profMean=NULL,typeProf="approx"){

  d<-objectUQ$kmModel@d

  if(is.null(plot_options)){
    plot_options<-list(save=F)
  }else{
    if(plot_options$save==TRUE && is.null(plot_options$folderPlots))
      plot_options$folderPlots <- './'

    if(is.null(plot_options$titleProf))
      plot_options$titleProf<-"Coordinate profiles"

    if(is.null(plot_options$title2d))
      plot_options$title2d<-"Posterior mean"

    if(is.null(plot_options$design)){
      plot_options$design<-matrix(NA,ncol=objectUQ$kmModel@d,nrow=100)
      for(i in seq(objectUQ$kmModel@d)){
        plot_options$design[,i]<-seq(0,1,,100)
      }
    }
    # Useful for serial plots in HPC
    if(is.null(plot_options$id_save)){
      plot_options$id_save<-""
    }

  }

  if(plot_options$save)
    cairo_pdf(filename = paste(plot_options$folderPlots,nameFile,plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
  mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
  par(mfrow=mfrows, mar = c(4, 5, 3, 1) + 0.1)
  for(coord in seq(d)){
    # create the title of each plot
    title_string<-paste("Coordinate",colnames(objectUQ$kmModel@X)[coord])
    # determine the plot lims
    if(is.null(plot_options$ylim)){
      if(typeProf=="approx"){
        ylimTemp<-range(c(range(objectUQ$profSups[coord,,]),range(objectUQ$profInfs[coord,,])))
      }else if(typeProf=="full"){
        ylimTemp<-range(c(range(objectUQ$profSups_full[coord,,]),range(objectUQ$profInfs_full[coord,,])))
      }else{
        ylimTemp<-range(c(range(objectUQ$profSups[coord,,]),range(objectUQ$profInfs[coord,,]),
                          range(objectUQ$profSups_full[coord,,]),range(objectUQ$profInfs_full[coord,,])))
      }
      if(!is.null(profMean)){
        ylimTemp<-range(c(ylimTemp,
                          range(profMean$res$min[,coord]),range(profMean$res$max[,coord])))
        #ylimTemp<-c(min(profMean$res$min[,coord]),max(profMean$res$max[,coord]))
        plot(plot_options$design[,coord],profMean$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      }else{
        ylimTemp<-range(c(ylimTemp,
                          range(objectUQ$prof_quantiles_approx$`0.5`$res$min[,coord]),
                          range(objectUQ$prof_quantiles_approx$`0.5`$res$max[,coord])))
        plot(plot_options$design[,coord],objectUQ$prof_quantiles_approx$`0.5`$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      }
    }else{
      ylimTemp<-plot_options$ylim[coord,]
    }
    for(i in seq(nsims)){
      if(typeProf=="approx"){
        dd<-seq(0,1,,length.out = length(objectUQ$profInfs[coord,,i]))
        lines(dd, objectUQ$profSups[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
        lines(dd, objectUQ$profInfs[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
      }else if (typeProf=="full"){
        dd<-seq(0,1,,length.out = length(objectUQ$profInfs_full[coord,,i]))
        lines(dd, objectUQ$profSups_full[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
        lines(dd, objectUQ$profInfs_full[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
      }else{
        dd<-seq(0,1,,length.out = length(objectUQ$profInfs[coord,,i]))
        lines(dd, objectUQ$profSups[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.3),lwd=0.2)
        lines(dd, objectUQ$profInfs[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.3),lwd=0.2)
        lines(dd, objectUQ$profSups_full[coord,,i],col=adjustcolor('darkorange3',alpha.f=0.3),lwd=0.2)
        lines(dd, objectUQ$profInfs_full[coord,,i],col=adjustcolor('darkorange3',alpha.f=0.3),lwd=0.2)
      }
    }
    if(plot_options$qq_fill){
      if(typeProf=="approx"){
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_approx[[2]]$res$min[,coord],
                  rev(objectUQ$prof_quantiles_approx[[1]]$res$min[,coord])),
                col = adjustcolor("red",alpha.f=0.5), border = NA)
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_approx[[2]]$res$max[,coord],
                  rev(objectUQ$prof_quantiles_approx[[1]]$res$max[,coord])),
                col = adjustcolor("red",alpha.f=0.5), border = NA,alpha=0.5)
      }else{
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_full[[2]]$res$min[,coord],
                  rev(objectUQ$prof_quantiles_full[[1]]$res$min[,coord])),
                col = adjustcolor("red",alpha.f=0.5), border = NA)
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_full[[2]]$res$max[,coord],
                  rev(objectUQ$prof_quantiles_full[[1]]$res$max[,coord])),
                col = adjustcolor("red",alpha.f=0.5), border = NA,alpha=0.5)
      }
    }
    #  lines(plot_options$design[,coord],allRes_approx$res$max[,coord],lty=3,col=4,lwd=2)
    #  points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$max[,coord],col=4)
    #  lines(plot_options$design[,coord],allRes_approx$res$min[,coord],lty=4,col=4,lwd=2)
    #  points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$min[,coord],col=4)
    ## @@@@
    ### NB changePP to be implemented!!!!!!!
    ## @@@
    #  abline(v=changePP$neverEx[[coord]],col=3,lwd=2.5)
    abline(h = threshold,col=2,lwd=2)
    if(typeProf=="approx"){
      for(j in seq(length(quantiles_uq))){
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$max[,coord],lty=j,col=1,lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$min[,coord],lty=j,col=1,lwd=1.5)
        #     abline(v=ccPP[[j]]$neverEx[[coord]],col=4,lwd=2,lty=2)
        #     abline(v=ccPP[[j]]$neverEx[[coord]],col=4,lwd=2,lty=3)
      }
    }else if(typeProf=="full"){
      if(!is.null(objectUQ$prof_quantiles_full)){
        for(j in seq(length(quantiles_uq))){
          lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$max[,coord],lty=j,col=adjustcolor('brown3',alpha.f=0.5),lwd=1.5)
          lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$min[,coord],lty=j,col=adjustcolor('brown3',alpha.f=0.5),lwd=1.5)
          #        abline(v=ccPP_full[[j]]$neverEx[[coord]],col="green3",lwd=2,lty=2)
          #        abline(v=ccPP_full[[j]]$neverEx[[coord]],col="green3",lwd=2,lty=3)
        }
      }
    }else{
      for(j in seq(length(quantiles_uq))){
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$max[,coord],lty=j,col=1,lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$min[,coord],lty=j,col=1,lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$max[,coord],lty=j,col=adjustcolor('brown3',alpha.f=0.5),lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$min[,coord],lty=j,col=adjustcolor('brown3',alpha.f=0.5),lwd=1.5)
        #     abline(v=ccPP[[j]]$neverEx[[coord]],col=4,lwd=2,lty=2)
        #     abline(v=ccPP[[j]]$neverEx[[coord]],col=4,lwd=2,lty=3)
      }
    }
    if(!is.null(profMean)){
      lines(plot_options$design[,coord],profMean$res$max[,coord],lty=1)
      lines(plot_options$design[,coord],profMean$res$min[,coord],lty=1)
    }else{
      lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx$`0.5`$res$max[,coord],lty=1)
      lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx$`0.5`$res$min[,coord],lty=1)
    }
    legend("bottomleft",c(as.expression(substitute(paste(P[coord]^sup,f,"/",P[coord]^inf,f," (full)"),list(coord=coord))),
                          as.expression(substitute(paste(P[coord]^sup,f[95]),list(coord=coord))),
                          as.expression(substitute(paste(P[coord]^inf,f[95]),list(coord=coord))),
                          as.expression(substitute(paste(P[coord]^sup,f[05]),list(coord=coord))),
                          as.expression(substitute(paste(P[coord]^inf,f[05]),list(coord=coord))),
                          "threshold"),
           lty=c(1,2,2,3,3,1),col=c(1,1,1,1,1,2),cex=1.2,lwd=c(2,2,2,2,2,2))
  }
  if(plot_options$save)
    dev.off()

}
