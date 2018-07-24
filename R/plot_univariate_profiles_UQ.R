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

  d <- ncol(plot_options$design)

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
      if(!is.null(objectUQ$bound$bound)){
        ylimTemp<-range(c(ylimTemp,
                          range(objectUQ$bound$bound$lower$res$min[,coord],na.rm = TRUE),
                          range(objectUQ$bound$bound$upper$res$min[,coord],na.rm = TRUE),
                          range(objectUQ$bound$bound$lower$res$max[,coord],na.rm = TRUE),
                          range(objectUQ$bound$bound$upper$res$max[,coord],na.rm = TRUE)))
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
 #       dd<-seq(0,1,,length.out = length(objectUQ$profInfs[coord,,i]))
        lines(plot_options$design[,coord], objectUQ$profSups[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
        lines(plot_options$design[,coord], objectUQ$profInfs[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
      }else if (typeProf=="full"){
 #       dd<-seq(0,1,,length.out = length(objectUQ$profInfs_full[coord,,i]))
        lines(plot_options$design[,coord], objectUQ$profSups_full[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
        lines(plot_options$design[,coord], objectUQ$profInfs_full[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
      }else{
  #      dd<-seq(0,1,,length.out = length(objectUQ$profInfs[coord,,i]))
        lines(plot_options$design[,coord], objectUQ$profSups[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.3),lwd=0.2)
        lines(plot_options$design[,coord], objectUQ$profInfs[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.3),lwd=0.2)
        lines(plot_options$design[,coord], objectUQ$profSups_full[coord,,i],col=adjustcolor('darkorange3',alpha.f=0.3),lwd=0.2)
        lines(plot_options$design[,coord], objectUQ$profInfs_full[coord,,i],col=adjustcolor('darkorange3',alpha.f=0.3),lwd=0.2)
      }
    }
    if(plot_options$qq_fill){
      if(typeProf=="approx"){
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_approx[[2]]$res$min[,coord],
                  rev(objectUQ$prof_quantiles_approx[[1]]$res$min[,coord])),
                col = plot_options$qq_fill_colors$approx, border = NA)
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_approx[[2]]$res$max[,coord],
                  rev(objectUQ$prof_quantiles_approx[[1]]$res$max[,coord])),
                col = plot_options$qq_fill_colors$approx, border = NA)
      }else{
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_full[[2]]$res$min[,coord],
                  rev(objectUQ$prof_quantiles_full[[1]]$res$min[,coord])),
                col = plot_options$qq_fill_colors$approx, border = NA)
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$prof_quantiles_full[[2]]$res$max[,coord],
                  rev(objectUQ$prof_quantiles_full[[1]]$res$max[,coord])),
                col = plot_options$qq_fill_colors$approx, border = NA)
      }
    }
    #  lines(plot_options$design[,coord],allRes_approx$res$max[,coord],lty=3,col=4,lwd=2)
    #  points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$max[,coord],col=4)
    #  lines(plot_options$design[,coord],allRes_approx$res$min[,coord],lty=4,col=4,lwd=2)
    #  points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$min[,coord],col=4)
    ## @@@@
    ### NB changePP to be implemented!!!!!!!
    if(!is.null(profMean)){
      changePP<-getChangePoints(threshold = threshold,allRes = profMean,Design = profMean$Design)
    }

    if(typeProf=="approx" || typeProf=="both"){
    ccPP_approx<-list()
    for(j in seq(length(quantiles_uq))){
      ccPP_approx[[j]]<-getChangePoints(threshold = threshold,allRes = objectUQ$prof_quantiles_approx[[j]],objectUQ$Design_approx)
    }
    names(ccPP_approx)<-quantiles_uq
    }

    if(typeProf=="full" || typeProf=="both"){
      ccPP_full<-list()
      for(j in seq(length(quantiles_uq))){
        ccPP_full[[j]]<-getChangePoints(threshold = threshold,allRes = objectUQ$prof_quantiles_full[[j]],objectUQ$Design_full)
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
      for(i in seq(num_T)){
        if(!is.null(profMean)){
          abline(v=changePP$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5)
          abline(v=changePP$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5)
        }else{
          abline(v=ccPP_approx[[j]]$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5,lty=1)
          abline(v=ccPP_approx[[j]]$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5,lty=1)
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
        for(i in seq(num_T)){
          if(!is.null(profMean)){
            abline(v=changePP$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5)
            abline(v=changePP$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5)
          }else{
            abline(v=ccPP_full[[j]]$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5)
            abline(v=ccPP_full[[j]]$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5)
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
      for(i in seq(num_T)){
        if(!is.null(profMean)){
          abline(v=changePP$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5)
          abline(v=changePP$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5)
        }else{
          abline(v=ccPP_full[[j]]$alwaysEx[[i]][[coord]],col=plot_options$col_CCPthresh_alw[i],lwd=2.5)
          abline(v=ccPP_full[[j]]$neverEx[[i]][[coord]],col=plot_options$col_CCPthresh_nev[i],lwd=2.5)
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
    if(plot_options$fun_evals>1){
      points(objectUQ$kmModel@X[,coord],objectUQ$kmModel@y,pch=17)
    }
    # bound
    if(!is.null(objectUQ$bound$bound)){
      lines(plot_options$design[,coord],objectUQ$bound$bound$lower$res$min[,coord],col=plot_options$bound_cols[1],lty=2,lwd=1.5)
      lines(plot_options$design[,coord],objectUQ$bound$bound$lower$res$max[,coord],col=plot_options$bound_cols[2],lty=2,lwd=1.5)
      lines(plot_options$design[,coord],objectUQ$bound$bound$upper$res$min[,coord],col=plot_options$bound_cols[1],lty=3,lwd=1.5)
      lines(plot_options$design[,coord],objectUQ$bound$bound$upper$res$max[,coord],col=plot_options$bound_cols[2],lty=3,lwd=1.5)

      if(plot_options$qq_fill){
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$bound$bound$upper$res$min[,coord],
                  rev(objectUQ$bound$bound$lower$res$min[,coord])),
                col = plot_options$qq_fill_colors$bound_min, border = NA)
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(objectUQ$bound$bound$upper$res$max[,coord],
                  rev(objectUQ$bound$bound$lower$res$max[,coord])),
                col = plot_options$qq_fill_colors$bound_max, border = NA)
      }else{
        lines(plot_options$design[,coord],objectUQ$bound$approx$lower$res$min[,coord],col=5,lty=2,lwd=1.2)
        lines(plot_options$design[,coord],objectUQ$bound$approx$lower$res$max[,coord],col=5,lty=2,lwd=1.2)
        lines(plot_options$design[,coord],objectUQ$bound$approx$upper$res$min[,coord],col=5,lty=3,lwd=1.2)
        lines(plot_options$design[,coord],objectUQ$bound$approx$upper$res$max[,coord],col=5,lty=3,lwd=1.2)
      }
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
  par(mfrow=c(1,1))
  if(plot_options$save)
    dev.off()

}


#' @name setPlotOptions
#' @title Set-up the plot options when NULL
#' @description  Function to set-up plot options for \link{plot_univariate_profiles_UQ}, \link{plotBivariateProfiles}, \link{coordinateProfiles}, \link{coordProf_UQ}, \link{obliqueProfiles} and \link{obliqueProf_UQ}.
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
#' \item{\code{qq_fill:}}{if TRUE it fills the region between the first 2 quantiles in \code{quantiles_uq} and between the upper and lower bound in \code{objectUQ$bound$bound}, if \code{NULL}, it is set as \code{FALSE}.}
#' \item{\code{bound_cols:}}{a vector of two strings containing the names of the colors for upper and lower bound plots.}
#' \item{\code{qq_fill_colors:}}{a list containing the colors for qq_fill: \code{approx} for 2 quantiles, \code{bound_min} for bounds on the profile inf, \code{bound_max} for profile sup. Initialized only if \code{qq_fill==TRUE}.}
#' \item{\code{col_CCPthresh_nev:}}{Color palette of dimension \code{num_T} for the colors of the vertical lines delimiting the intersections between the profiles sup and the thresholds}
#' \item{\code{col_CCPthresh_alw:}}{Color palette of dimension \code{num_T} for the colors of the vertical lines delimiting the intersections between the profiles inf and the thresholds}
#' \item{\code{col_thresh:}}{Color palette of dimension \code{num_T} for the colors of the thresholds}
#' \item{\code{fun_evals:}{integer denoting the level of plot for the true evaluations. \itemize{
#' \item{0: }{default, no plots for true evaluations;}
#' \item{1: }{plot the true evaluations as points in 2d plots, no true evaluation plots in 1d};
#' \item{2: }{plot true evaluations, in 2d with different color for values above threshold;}
#' \item{3: }{plot true evaluations, in 2d plots in color, with background of the image colored as proportion of points inside excursion;} }}}
#' }
#' if all the fields are already filled then returns \code{plot_options}
#' @export
setPlotOptions<-function(plot_options=NULL,d,num_T,kmModel=NULL){

  # Start with saving, folder options, titles
  full_plot_options <-list(save=F, folderPlots = './', id_save="",
                           titleProf="Profile extrema", title2d="Posterior mean")



  # Design plot
  full_plot_options$design<-matrix(NA,ncol=d,nrow=100)
  for(i in seq(d)){
    full_plot_options$design[,i]<-seq(0,1,,100)
  }


  # coordinate names
  if(is.null(kmModel)){
    full_plot_options$coord_names<-apply(matrix(1:d,ncol=1),1,function(x) paste("x_",x,sep=''))
  }else{
    full_plot_options$coord_names<-colnames(kmModel@X)
  }

  ## Colors set-up

  # AlwaysEx colors
  full_plot_options$col_CCPthresh_alw<-adjustcolor(RColorBrewer::brewer.pal(n=max(num_T+1,3),name='Purples'),offset = c(-0.25, -0.25, -0.25, 0))
  full_plot_options$col_CCPthresh_alw<-full_plot_options$col_CCPthresh_alw[(2:length(full_plot_options$col_CCPthresh_alw))]
  if(num_T==1)
    full_plot_options$col_CCPthresh_alw[2]

  # NeverEx colors
  full_plot_options$col_CCPthresh_nev<-adjustcolor(RColorBrewer::brewer.pal(n=max(num_T+1,3),name='Greens'),offset = c(-0.25, -0.25, -0.25, 0))
  full_plot_options$col_CCPthresh_nev<-full_plot_options$col_CCPthresh_nev[(2:length(full_plot_options$col_CCPthresh_nev))]
  if(num_T==1)
    full_plot_options$col_CCPthresh_nev[2]

  # threshold colors
  full_plot_options$col_thresh<-adjustcolor(RColorBrewer::brewer.pal(n=max(num_T+1,3),name='Reds'),offset = c(-0.3, -0.3, -0.3, 0))
  full_plot_options$col_thresh<-full_plot_options$col_thresh[(2:length(full_plot_options$col_thresh))]
  if(num_T==1)
    full_plot_options$col_thresh[2]

  # bound colors
  # first one is min, second one is max
  full_plot_options$bound_cols<-c("darkseagreen4","deepskyblue4")



  # plot function evaluations
  full_plot_options$fun_evals <- 0

  # qq_fill
  full_plot_options$qq_fill <- FALSE
  # qq_fill colors
  full_plot_options$qq_fill_colors<-list(
    approx=adjustcolor("red",alpha.f=0.5),
    bound_max=adjustcolor(full_plot_options$bound_cols[2],alpha.f=0.2),
    bound_min=adjustcolor(full_plot_options$bound_cols[1],alpha.f=0.2) )

  if(is.null(plot_options)){
    plot_options<-full_plot_options
  }else{
    plot_options <- modifyList(full_plot_options,plot_options)
  }
  return(plot_options)
}
