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

  num_T<-length(threshold)

  plot_options<-setPlotOptions(plot_options = plot_options,d=d,num_T=num_T,kmModel = objectUQ$kmModel)


  if(plot_options$save)
    cairo_pdf(filename = paste(plot_options$folderPlots,nameFile,plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
  mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
  par(mfrow=mfrows, mar = c(4, 5, 3, 1) + 0.1)
  for(coord in seq(d)){
    # create the title of each plot
    title_string<-paste("Coordinate",plot_options$coord_names[coord])
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
      if(!is.null(profMean)){
        plot(plot_options$design[,coord],profMean$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      }else{
        plot(plot_options$design[,coord],objectUQ$prof_quantiles_approx$`0.5`$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      }
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
    if(!is.null(profMean)){
      changePP<-getChangePoints(threshold = threshold,allRes = profMean)
      for(i in seq(num_T)){
        abline(v=changePP$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5)
        abline(v=changePP$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5)
      }
    }

    if(typeProf=="approx" || typeProf=="both"){
    ccPP_approx<-list()
    for(j in seq(length(quantiles_uq))){
      ccPP_approx[[j]]<-getChangePoints(threshold = threshold,allRes = objectUQ$prof_quantiles_approx[[j]])
    }
    names(ccPP_approx)<-quantiles_uq
    }

    if(typeProf=="full" || typeProf=="both"){
      ccPP_full<-list()
      for(j in seq(length(quantiles_uq))){
        ccPP_full[[j]]<-getChangePoints(threshold = threshold,allRes = objectUQ$prof_quantiles_full[[j]])
      }
      names(ccPP_full)<-quantiles_uq
    }

    ## @@@
    abline(h = threshold,col=plot_options$col_thresh,lwd=2)
    if(typeProf=="approx"){
      for(j in seq(length(quantiles_uq))){
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$max[,coord],lty=j+1,col=1,lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$min[,coord],lty=j+1,col=1,lwd=1.5)
        for(tt in seq(num_T)){
          abline(v=ccPP_approx[[j]]$alwaysEx[[tt]][[coord]],col=plot_options$col_CCPthresh_alw[tt],lwd=2,lty=2)
          abline(v=ccPP_approx[[j]]$neverEx[[tt]][[coord]],col=plot_options$col_CCPthresh_nev[tt],lwd=2,lty=3)
        }
      }
    }else if(typeProf=="full"){
      if(!is.null(objectUQ$prof_quantiles_full)){
        for(j in seq(length(quantiles_uq))){
          lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$max[,coord],lty=j+1,col=1,lwd=1.5)
          lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$min[,coord],lty=j+1,col=1,lwd=1.5)
          for(tt in seq(num_T)){
            abline(v=ccPP_full[[j]]$alwaysEx[[tt]][[coord]],col=plot_options$col_CCPthresh_alw[tt],lwd=2,lty=2)
            abline(v=ccPP_full[[j]]$neverEx[[tt]][[coord]],col=plot_options$col_CCPthresh_nev[tt],lwd=2,lty=3)
          }
        }
      }
    }else{
      for(j in seq(length(quantiles_uq))){
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$max[,coord],lty=j+1,col=1,lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx[[j]]$res$min[,coord],lty=j+1,col=1,lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$max[,coord],lty=j+1,col=adjustcolor('brown3',alpha.f=0.5),lwd=1.5)
        lines(plot_options$design[,coord],objectUQ$prof_quantiles_full[[j]]$res$min[,coord],lty=j+1,col=adjustcolor('brown3',alpha.f=0.5),lwd=1.5)
        for(tt in seq(num_T)){
          abline(v=ccPP_approx[[j]]$alwaysEx[[tt]][[coord]],col=4,lwd=2,lty=2)
          abline(v=ccPP_approx[[j]]$neverEx[[tt]][[coord]],col=4,lwd=2,lty=3)
          abline(v=ccPP_full[[j]]$alwaysEx[[tt]][[coord]],col=adjustcolor('brown4',alpha.f=0.5),lwd=2,lty=3)
          abline(v=ccPP_full[[j]]$neverEx[[tt]][[coord]],col=adjustcolor('brown3',alpha.f=0.5),lwd=2,lty=3)
        }
      }
    }
    if(!is.null(profMean)){
      lines(plot_options$design[,coord],profMean$res$max[,coord],lty=1)
      lines(plot_options$design[,coord],profMean$res$min[,coord],lty=1)
    }else{
      lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx$`0.5`$res$max[,coord],lty=1)
      lines(plot_options$design[,coord],objectUQ$prof_quantiles_approx$`0.5`$res$min[,coord],lty=1)
    }
    if(plot_options$fun_evals){
      points(objectUQ$kmModel@X[,coord],objectUQ$kmModel@y,pch=17)
    }

    if(!is.null(plot_options$legend)){
      legend("bottomleft",c(as.expression(substitute(paste(P[coord]^sup,f,"/",P[coord]^inf,f," (full)"),list(coord=coord))),
                            as.expression(substitute(paste(P[coord]^sup,f[qq]),list(coord=coord,qq=quantiles_uq[1]))),
                            as.expression(substitute(paste(P[coord]^inf,f[qq]),list(coord=coord,qq=quantiles_uq[1]))),
                            as.expression(substitute(paste(P[coord]^sup,f[qq]),list(coord=coord,qq=quantiles_uq[2]))),
                            as.expression(substitute(paste(P[coord]^inf,f[qq]),list(coord=coord,qq=quantiles_uq[2]))),
                            "threshold"),
             lty=c(1,2,2,3,3,1),col=c(1,1,1,1,1,2),cex=1.2,lwd=c(2,2,2,2,2,2))
    }
  }
  if(plot_options$save)
    dev.off()

}


#' @name setPlotOptions
#' @title Set-up the plot options when NULL
#' @description  function to set-up plot options for plot_univariate_profiles_UQ, coordProf_UQ and coordinateProfiles
#' @author Dario Azzimonti
#' @param plot_options the list of plot options to set-up
#' @param d number of coordinates
#' @param num_T number of thresholds of interest
#' @param kmModel a \link[DiceKriging]{km} model, used to obtain the coordinates names.
#' @return the properly set-up list containing the following fields \itemize{
#' \item{\code{save:}}{boolean, if TRUE saves the plots in \code{folderPlots}}
#' \item{\code{folderPlots:}}{a string containing the destination folder for plots, if \code{save==TRUE} default is \code{./}}
#' \item{\code{ylim:}}{a matrix \code{coord}x2 containing the ylim for each coordinate, if NULL in \code{plot_options} this is left NULL and automatically set at the plot time.}
#' \item{\code{titleProf:}}{a string containing the title for the coordinate profile plots, default is \code{"Coordinate profiles"}}
#' \item{\code{title2d:}}{a string containing the title for the 2d plots (if the input is 2d), default is \code{"Posterior mean"}}
#' \item{\code{design:}}{a \eqn{dxr} matrix where \eqn{d} is the input dimension and \eqn{r} is the size of the discretization for plots at each dimension}
#' \item{\code{coord_names:}}{a \eqn{d}-vector of characters naming the dimensions. If NULL and \code{kmModel} not NULL then it is the names of \code{kmModel@X} otherwise \code{x_1,...,x_d}}
#' \item{\code{id_save:}}{a string to be added to the plot file names, useful for serial computations on HPC, left as in \code{plot_options}.}
#' \item{\code{qq_fill:}}{if TRUE it fills the region between the first 2 quantiles in \code{quantiles_uq}, left as in \code{plot_options}.}
#' \item{\code{col_CCPthresh_nev:}}{Color palette of dimension \code{num_T} for the colors of the vertical lines delimiting the intersections between the profiles sup and the thresholds}
#' \item{\code{col_CCPthresh_alw:}}{Color palette of dimension \code{num_T} for the colors of the vertical lines delimiting the intersections between the profiles inf and the thresholds}
#' \item{\code{col_thresh:}}{Color palette of dimension \code{num_T} for the colors of the thresholds}
#' }
#' if all the fields are already filled then returns \code{plot_options}
#' @export
setPlotOptions<-function(plot_options=NULL,d,num_T,kmModel=NULL){

  # Start with saving and folder options
  if(is.null(plot_options)){
    # don't save to file
    plot_options<-list(save=F)
  }else{
    # set save directory
    if(plot_options$save==TRUE && is.null(plot_options$folderPlots))
      plot_options$folderPlots <- './'

    # Useful for serial plots in HPC
    if(is.null(plot_options$id_save))
      plot_options$id_save<-""
  }

  ## Titles of plots
  # profile plots
  if(is.null(plot_options$titleProf))
    plot_options$titleProf<-"Coordinate profiles"
  # 2d plot
  if(is.null(plot_options$title2d))
    plot_options$title2d<-"Posterior mean"

  # Design plot
  if(is.null(plot_options$design)){
    plot_options$design<-matrix(NA,ncol=d,nrow=100)
    for(i in seq(d)){
      plot_options$design[,i]<-seq(0,1,,100)
    }
  }

  # coordinate names
  if(is.null(plot_options$coord_names)){
    if(is.null(kmModel)){
      plot_options$coord_names<-apply(matrix(1:d,ncol=1),1,function(x) paste("x_",x,sep=''))
    }else{
      plot_options$coord_names<-colnames(kmModel@X)
    }
  }

  ## Colors set-up

  # AlwaysEx colors
  if(is.null(plot_options$col_CCPthresh_alw)){
    plot_options$col_CCPthresh_alw<-adjustcolor(RColorBrewer::brewer.pal(n=max(num_T+1,3),name='Purples'),offset = c(-0.2, -0.2, -0.2, 0))
    plot_options$col_CCPthresh_alw<-plot_options$col_CCPthresh_alw[(2:length(plot_options$col_CCPthresh_alw))]
    if(num_T==1)
      plot_options$col_CCPthresh_alw[2]
  }
  # NeverEx colors
  if(is.null(plot_options$col_CCPthresh_nev)){
    plot_options$col_CCPthresh_nev<-adjustcolor(RColorBrewer::brewer.pal(n=max(num_T+1,3),name='Greens'),offset = c(-0.2, -0.2, -0.2, 0))
    plot_options$col_CCPthresh_nev<-plot_options$col_CCPthresh_nev[(2:length(plot_options$col_CCPthresh_nev))]
    if(num_T==1)
      plot_options$col_CCPthresh_nev[2]
  }
  # threshold colors
  if(is.null(plot_options$col_thresh)){
    plot_options$col_thresh<-adjustcolor(RColorBrewer::brewer.pal(n=max(num_T+1,3),name='Reds'),offset = c(-0.2, -0.2, -0.2, 0))
    plot_options$col_thresh<-plot_options$col_thresh[(2:length(plot_options$col_thresh))]
    if(num_T==1)
      plot_options$col_thresh[2]
  }

  # plot function evaluations
  if(is.null(plot_options$fun_evals))
    plot_options$fun_evals <- FALSE


  return(plot_options)
}
