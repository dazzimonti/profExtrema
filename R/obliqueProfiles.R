# obliqueProfiles function
#' @author Dario Azzimonti
#' @name obliqueProfiles
#' @title Oblique coordinate profiles starting from a kriging model
#' @description The function obliqueProfiles computes the (oblique) profile extrema functions for the posterior mean of a Gaussian process and its confidence bounds
#' @param object either a \link[DiceKriging]{km} model or a list containing partial results. If \code{object} is a km model then all computations are carried out. If \code{object} is a list, then the function carries out all computations to complete the list results.
#' @param allPsi a list containing the matrices Psi (dim \eqn{pxd}) for which to compute the profile extrema
#' @param threshold the threshold of interest
#' @param options_full an optional list of options for getProfileExtrema, see \link{getProfileExtrema} for details.
#' @param options_approx an optional list of options for approxProfileExtrema, see \link{approxProfileExtrema} for details.
#' @param uq_computations boolean, if TRUE the uq computations for the profile mean are computed.
#' @param plot_level an integer to select the plots to return (0=no plots, 1=basic plots, 2= all plots)
#' @param plot_options an optional list of parameters for plots. See \link{setPlotOptions} for currently available options.
#' @param CI_const an optional vector containing the constants for the CI. If not NULL, then profiles extrema for \eqn{m_n(x) \pm CI_const[i]*s_n(x,x)} are computed.
#' @param return_level an integer to select the amount of details returned
#' @param ... additional parameters to be passed to \link{obliqueProf_UQ}.
#' @return If return_level=1 a list containing
#' \itemize{
#' \item{\code{profMean_full:}}{the results of \code{getProfileExtrema} for the posterior mean}
#' \item{\code{profMean_approx:}}{the results of \code{approxProfileExtrema} for the posterior mean}
#' \item{\code{res_UQ:}}{the results of \code{obliqueProf_UQ} for the posterior mean}
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
#' @examples
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # Compute a kriging model from 50 evaluations of the Branin function
#' # Define the function
#' g=function(x){
#'   return(-branin(x))
#' }
#' gp_des<-lhs::maximinLHS(20,2)
#' reals<-apply(gp_des,1,g)
#' kmModel<-km(design = gp_des,response = reals,covtype = "matern3_2")
#'
#' threshold=-10
#'
#' # Compute oblique profiles on the posterior mean
#' # (for theta=0 it is equal to coordinateProfiles)
#' options_full<-list(multistart=4,heavyReturn=TRUE,discretization=100)
#' options_approx<- list(multistart=4,heavyReturn=TRUE,initDesign=NULL,fullDesignSize=100)
#' theta=pi/4
#' allPsi = list(Psi1=matrix(c(cos(theta),sin(theta)),ncol=2),
#' Psi2=matrix(c(cos(theta+pi/2),sin(theta+pi/2)),ncol=2))
#'
#' \dontrun{
#' profMeans<-obliqueProfiles(object = kmModel,allPsi = allPsi,threshold = threshold,
#'                            options_full = options_full,options_approx = options_approx,
#'                            uq_computations = FALSE,plot_level = 3,plot_options = NULL,
#'                            CI_const = NULL,return_level = 2)
#'
#'
#' # Approximate oblique profiles with UQ
#' plot_options<-list(save=FALSE, titleProf = "Coordinate profiles",
#'                    title2d = "Posterior mean",qq_fill=TRUE)
#' options_sims<-list(nsim=150)
#' obProfUQ<-obliqueProfiles(object=profMeans,threshold=threshold,allPsi = allPsi,
#'                            options_full=options_full, options_approx=options_approx,
#'                            uq_computations=TRUE, plot_level=3,plot_options=NULL,
#'                            CI_const=NULL,return_level=2,options_sims=options_sims)
#' }
#' @export
obliqueProfiles = function(object,allPsi,threshold,options_full=NULL,options_approx=NULL,uq_computations=FALSE,plot_level=0,plot_options=NULL,CI_const=NULL,return_level=1,...){

  # Check object
  if(is(object,"km")){
    object<-list(kmModel=object)
  }else if(!is.list(object)){
    stop("object must be either a list or a km object")
  }

  # save the dimension of the input
  d<-object$kmModel@d

  # Options setup
  if(is.null(options_full)){
    options_full<-list(multistart=8,heavyReturn=TRUE,discretization=100)
  }

  # Initialize options_approx
  if(is.null(options_approx)){
    options_approx<- list(multistart=8,heavyReturn=TRUE,fullDesignSize=100,smoother=NULL)
  }

  # initialize limits
  if(is.null(options_approx$limits)){
    ll_b = rep(0,d)
    uu_b = rep(1,d)
  }else{
    ll_b = options_approx$limits$lower
    uu_b = options_approx$limits$upper
  }
  # Useful for choosing limits for eta
  cubeVertex<-matrix(c(ll_b[1],uu_b[1]),nrow=1)
  ii=1
  while(ii<d){
    ii=ii+1
    cubeVertex <- cbind(cubeVertex,cubeVertex)
    cubeVertex <-rbind(cubeVertex, c(rep(ll_b[ii],ncol(cubeVertex)/2),c(rep(uu_b[ii],ncol(cubeVertex)/2))))
  }

  # Initialize initial design for approximation
  num_Psi = length(allPsi)
  if(is.null(options_approx$initDesign)){
    init_des<-list()
    for(i in seq(num_Psi)){
      p<-nrow(allPsi[[i]])
      # Choose limits for etas for current Psi
      mmEtas<-apply(crossprod(t(allPsi[[i]]),cubeVertex),1,min)
      MMetas<-apply(crossprod(t(allPsi[[i]]),cubeVertex),1,max)

      # Get initial design
      if(length(options_approx$initDesign)<num_Psi){
        if(p==1){
          init_des[[i]]<-matrix(seq(from=mmEtas[1],to=MMetas[1],,ceiling(sqrt(d)*10)),ncol=1)
        }else{
          ### NEEDS TEST!
          init_des[[i]]<-matrix(mmEtas,ncol=p,byrow = T,nrow=ceiling(sqrt(d*p)*10))+maximinLHS(ceiling(sqrt(d*p)*10),p)%*%diag(MMetas-mmEtas,ncol=p)
        }
      }
  #    init_des<-maximinLHS(ceiling(sqrt(object$kmModel@d)*10),object$kmModel@d)
    }
    options_approx$initDesign<-init_des
  }
  p<-nrow(allPsi[[1]])

  # save number of thresholds
  num_T<-length(threshold)

  # Set up plot options
  plot_options<-setPlotOptions(plot_options = plot_options,d=d,num_T=num_T,kmModel=object$kmModel)

  ## posterior mean part ##
  # Let us define the posterior mean function and gradient in a 'optimizer-friendly' way
  g <-function(y){
    y<-matrix(y,ncol=d)
    colnames(y)<-colnames(object$kmModel@X)
    return(predict.km(object = object$kmModel,newdata = y,type="UK",light.return = TRUE,se.compute = FALSE)$mean)
  }
  gprime <- function(y){
    y<-matrix(y,ncol=d)
    colnames(y)<-colnames(object$kmModel@X)
    return(t(gradKm_dnewdata(object = object$kmModel,newdata = y,type = "UK",light.return = TRUE,se.compute = FALSE)$mean))
  }

  # Compute profile mean full, if not already computed
  if(is.null(object$profMean_full)){
    timeIn<-get_nanotime()
    object$profMean_full<-getProfileExtrema(f = g,fprime = gprime,d = d,allPsi = allPsi,opts = options_full)
    timeMean_full<-(get_nanotime()-timeIn)*1e-9
  }else{
    if(return_level>1 && is.null(object$more$times$full))
      object$more<-list(times=list(full=NA))
  }
  #    profMean_full,profMean_approx,res_UQ
  #    times<-list(full=timeMean_full,approx=timeMean_approx)
  #    results$more<-list(abs_err=abs_err,times=times)
  # Compute profile extrema with full optim

  gc()

  # Obtain the threshold points selected by profile extrema
  if(p==1){
    changePP<-getChangePoints(threshold = threshold,allRes = object$profMean_full,Design = object$profMean_full$Design)
  }else{
    changePP<-NULL
  }

  # Plot posterior mean, ONLY for dimension==2
  if(plot_level>=2 & num_Psi==2 & d==2){

    # since dimension==2 we can plot the posterior mean
    newdata<-expand.grid(seq(ll_b[1],uu_b[1],,100),seq(ll_b[2],uu_b[2],,100))
    colnames(newdata)<-colnames(object$kmModel@X)
    pred2d<-predict.km(object$kmModel,newdata = newdata,type = "UK",light.return = TRUE,se.compute = FALSE)

    # Plot the posterior mean and visualize the actual excursion set and the regions of no-excursion according to the profile extrema functions.
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"2dpostMean_1",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    par(mar = c(5, 5, 4, 2) + 0.1,mfrow=c(1,1))
    image(matrix(pred2d$mean,nrow = 100),col=gray.colors(20), main=plot_options$title2d,xlab = plot_options$coord_names[1],ylab= plot_options$coord_names[2],
          cex.main=3,cex.axis=1.8,cex.lab=2.8)
    contour(matrix(pred2d$mean,nrow = 100),add=T,nlevels = 10,lwd=1.5,labcex=1.2)
    contour(matrix(pred2d$mean,nrow = 100),add=T,levels = threshold,col=plot_options$col_thresh,lwd=3,labcex=1.5)
    for(tt in seq(num_T)){
      plotOblique(changePoints = changePP$alwaysEx[[tt]][[1]],direction = allPsi[[1]],col=plot_options$col_CCPthresh_alw[tt],lwd=2.5)
      plotOblique(changePoints = changePP$alwaysEx[[tt]][[2]],direction = allPsi[[2]],col=plot_options$col_CCPthresh_alw[tt],lwd=2.5)
      plotOblique(changePoints = changePP$neverEx[[tt]][[1]],direction = allPsi[[1]],col=plot_options$col_CCPthresh_nev[tt],lwd=2.5)
      plotOblique(changePoints = changePP$neverEx[[tt]][[2]],direction = allPsi[[2]],col=plot_options$col_CCPthresh_nev[tt],lwd=2.5)
    }
    if(plot_options$fun_evals>0){
      points(object$kmModel@X,pch=17,cex=1.6)
    }
    if(plot_options$save)
      dev.off()
  }

  # Plot the profile extrema functions
  if(plot_level>0){
    if(p==1){
      if(plot_options$save)
        cairo_pdf(filename = paste(plot_options$folderPlots,"profMean_full",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
      colnames(object$profMean_full$res$min)<-plot_options$coord_names
      plotMaxMin(allRes = object$profMean_full,threshold = threshold,Design = object$profMean_full$Design)
    }else{
      if(plot_options$save)
        pdf(file = paste(plot_options$folderPlots,"profMean_full",plot_options$id_save,".pdf",sep=""),width = 18,height = 9)
      plotOneBivProfile(allRes = object$profMean_full,allPsi = allPsi,Design=object$profMean_full$Design,threshold=threshold,main_addendum="(mean, full)",xlab=expression(eta[1]),ylab=expression(eta[2]))
    }
    if(plot_options$save)
      dev.off()
  }

  # Compute profile extrema with approximation, if not already computed
  if(is.null(object$profMean_approx)){
    timeIn<-get_nanotime()
    object$profMean_approx<-approxProfileExtrema(f = g,fprime = gprime,d = d,allPsi = allPsi,opts = options_approx)
    timeMean_approx<-(get_nanotime()-timeIn)*1e-9
  }else{
    if(return_level>1 && is.null(object$more$times$approx))
      object$more$times$approx=NA
  }

  if(plot_level>0 & p==1){
    oldpar<-par()
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"profMean_comparison",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
    par(mfrow=mfrows, mar = c(4, 5, 3, 1) + 0.1)
    for(coord in seq(num_Psi)){
      if(is.null(plot_options$ylim)){
        ylimTemp<-c(min(object$profMean_full$res$min[,coord]),max(object$profMean_full$res$max[,coord]))
      }else{
        ylimTemp<-plot_options$ylim[coord,]
      }
      # create the title of each plot
      title_string<-paste("Coordinate",plot_options$coord_names[coord])

      plot(object$profMean_full$Design[,coord],object$profMean_full$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
           xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      lines(object$profMean_full$Design[,coord],object$profMean_full$res$max[,coord],lty=1)
      lines(object$profMean_approx$Design[,coord],object$profMean_approx$res$max[,coord],lty=3,col=4,lwd=2)
      points(t(object$profMean_approx$profPoints$design[[coord]]),object$profMean_approx$profPoints$e_res$max[,coord],col=4)
      lines(object$profMean_approx$Design[,coord],object$profMean_approx$res$min[,coord],lty=4,col=4,lwd=2)
      points(t(object$profMean_approx$profPoints$design[[coord]]),object$profMean_approx$profPoints$e_res$min[,coord],col=4)
      for(tt in seq(num_T)){
        abline(v=changePP$neverEx[[tt]][[coord]],col=plot_options$col_CCPthresh_nev[tt],lwd=2.5)
        abline(v=changePP$alwaysEx[[tt]][[coord]],col=plot_options$col_CCPthresh_alw[tt],lwd=2.5)
      }
      abline(h = threshold,col=plot_options$col_thresh,lwd=2)
      if(!is.null(plot_options$legend)){
        legend("bottomleft",c(as.expression(substitute(paste(P[coord]^sup,f," (full)"),list(coord=coord))),
                              as.expression(substitute(paste(P[coord]^inf,f," (full)"),list(coord=coord))),
                              as.expression(substitute(paste(P[coord]^sup,f," (1order)"),list(coord=coord))),
                              as.expression(substitute(paste(P[coord]^inf,f," (1order)"),list(coord=coord))),
                              "threshold"),
               lty=c(1,1,3,4,1),col=c(1,1,4,4,2),cex=1.8,lwd=c(2,2,2,2,2))
      }
    }
    par(mfrow=c(1,1))
    if(plot_options$save)
      dev.off()
    #  par(oldpar)
  }else{
    if(plot_options$save)
      pdf(file = paste(plot_options$folderPlots,"profMean_approx",plot_options$id_save,".pdf",sep=""),width = 18,height = 9)
    plotOneBivProfile(allRes = object$profMean_approx,allPsi = allPsi,Design=object$profMean_approx$Design,threshold=threshold,main_addendum="(mean, approx)",xlab=expression(eta[1]),ylab=expression(eta[2]))
    if(plot_options$save)
      dev.off()
  }

  # save absolute errors with approximation
  if(nrow(object$profMean_approx$res$min)==nrow(object$profMean_full$res$min)){
  abs_err_min<-max(abs(object$profMean_approx$res$min-object$profMean_full$res$min)/object$profMean_full$res$min)
  abs_err_max<-max(abs(object$profMean_approx$res$max-object$profMean_full$res$max)/object$profMean_full$res$max)
  abs_err<-c(abs_err_min,abs_err_max)
  }else{
    warning("Different discretizations for approx and full profiles, error not computed!")
    abs_err = c(NA,NA)
  }
#   abs_err = c(NA,NA)
  ## End posterior mean part ##

  ## quantiles part ##
  if(!is.null(CI_const)){
    g_s2<-function(y){
      y<-matrix(y,ncol=d)
      colnames(y)<-colnames(object$kmModel@X)
      return(predict.km(object = object$kmModel,newdata = y,type="UK",light.return = TRUE,se.compute = TRUE)$sd)
    }
    g_s2prime <- function(y){
      y<-matrix(y,ncol=d)
      colnames(y)<-colnames(object$kmModel@X)
      s2<-predict.km(object = object$kmModel,newdata = y,type="UK",light.return = TRUE,se.compute = TRUE)$sd
      s2grad<-gradKm_dnewdata(object = object$kmModel,newdata = y,type = "UK",se.compute = TRUE)$s2
      return(0.5*t(s2grad)/s2)
    }

    # Compute profile sd full, if not already computed
    if(is.null(object$profSd_full)){
      timeIn<-get_nanotime()
      object$profSd_full<-getProfileExtrema(f = g_s2,fprime = g_s2prime,d = d,allPsi = allPsi,opts = options_full)
      timeVar_full<-(get_nanotime()-timeIn)*1e-9
    }else{
      if(return_level>1 && is.null(object$more$times$full))
        object$more$times$full_var=NA
    }

    # Obtain the threshold points selected by profile extrema
    prof_CI_const_full<-list()
    temp<-list()
    changePP_CI_full<-list()
    for(i in seq(length(CI_const))){
      temp_m<-list(res=list())
      temp_m$res<-list(min= object$profMean_full$res$min - CI_const[i]*object$profSd_full$res$min,
                       max= object$profMean_full$res$max - CI_const[i]*object$profSd_full$res$max)
      temp_M<-list(res=list())
      temp_M$res<-list(min= object$profMean_full$res$min + CI_const[i]*object$profSd_full$res$min,
                       max= object$profMean_full$res$max + CI_const[i]*object$profSd_full$res$max)
      prof_CI_const_full[[i]]<-list(lower=temp_m, upper=temp_M)

      if(p==1)
        changePP_CI_full[[i]]<-list(lower=getChangePoints(threshold = threshold,allRes = temp_m,Design = object$profMean_full$Design),
                                  upper=getChangePoints(threshold = threshold,allRes = temp_M,Design = object$profMean_full$Design))
    }

    # Compute profile extrema var with approximation, if not already computed
    #    if(is.null(object$profVar_approx)){
    #      timeIn<-get_nanotime()
    #      object$profVar_approx<-approxMaxMin(f = g_s2,fprime = g_s2prime,d = d,opts = options_approx)
    #      timeVar_approx<-(get_nanotime()-timeIn)*1e-9
    #    }else{
    #      if(return_level>1 && is.null(object$more$times$approx))
    #        object$more$times$approx=NA
    #    }

    if(plot_level>0 & p==1){
      oldpar<-par()
      if(plot_options$save)
        cairo_pdf(filename = paste(plot_options$folderPlots,"profMeanCI_comparison",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
      mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
      par(mfrow=mfrows, mar = c(4, 5, 3, 1) + 0.1)
      for(coord in seq(num_Psi)){
        if(is.null(plot_options$ylim)){
          ylimTemp_CI<-c(object$profMean_full$res$min[,coord],object$profMean_full$res$max[,coord])
          for(i in seq(length(CI_const))){
            ylimTemp_CI<-c(ylimTemp_CI,prof_CI_const_full[[i]]$lower$res$min[,coord],
                           prof_CI_const_full[[i]]$lower$res$max[,coord],
                           prof_CI_const_full[[i]]$upper$res$min[,coord],
                           prof_CI_const_full[[i]]$upper$res$max[,coord])
          }
          ylimTemp_CI<-range(ylimTemp_CI)
        }else{
          ylimTemp<-plot_options$ylim[coord,]
        }

        # create the title of each plot
        title_string<-paste("Coordinate",plot_options$coord_names[coord])
        # full mean
        plot(object$profMean_full$Design[,coord],object$profMean_full$res$min[,coord],ylim=ylimTemp_CI,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
        lines(object$profMean_full$Design[,coord],object$profMean_full$res$max[,coord],lty=1)
        # approx mean
        lines(object$profMean_approx$Design[,coord],object$profMean_approx$res$max[,coord],lty=3,col=4,lwd=2)
        points(t(object$profMean_approx$profPoints$design[[coord]]),object$profMean_approx$profPoints$e_res$max[,coord],col=4)
        lines(object$profMean_approx$Design[,coord],object$profMean_approx$res$min[,coord],lty=4,col=4,lwd=2)
        points(t(object$profMean_approx$profPoints$design[[coord]]),object$profMean_approx$profPoints$e_res$min[,coord],col=4)
        for(tt in seq(num_T)){
          abline(v=changePP$neverEx[[tt]][[coord]],col=plot_options$col_CCPthresh_nev[tt],lwd=2.5)
          abline(v=changePP$alwaysEx[[tt]][[coord]],col=plot_options$col_CCPthresh_alw[tt],lwd=2.5)
        }
        # full CI_const
        CI_cols1<-heat.colors(n = length(CI_const)+1,alpha=0.7)
        CI_cols2<-terrain.colors(n = length(CI_const)+1,alpha=0.7)
        for(i in seq(length(CI_const))){
          lines(object$profMean_full$Design[,coord],prof_CI_const_full[[i]]$upper$res$max[,coord],lty=4,col=CI_cols1[i])
          lines(object$profMean_full$Design[,coord],prof_CI_const_full[[i]]$lower$res$max[,coord],lty=4,col=CI_cols1[i])
          lines(object$profMean_full$Design[,coord],prof_CI_const_full[[i]]$upper$res$min[,coord],lty=4,col=CI_cols2[i])
          lines(object$profMean_full$Design[,coord],prof_CI_const_full[[i]]$lower$res$min[,coord],lty=4,col=CI_cols2[i])
        }
        abline(h = threshold,col=plot_options$col_thresh,lwd=2)
        if(!is.null(plot_options$legend)){
          legend("bottomleft",c(as.expression(substitute(paste(P[coord]^sup,f," (full)"),list(coord=coord))),
                                as.expression(substitute(paste(P[coord]^inf,f," (full)"),list(coord=coord))),
                                as.expression(substitute(paste(P[coord]^sup,f," (1order)"),list(coord=coord))),
                                as.expression(substitute(paste(P[coord]^inf,f," (1order)"),list(coord=coord))),
                                "threshold"),
                 lty=c(1,1,3,4,1),col=c(1,1,4,4,2),cex=1.8,lwd=c(2,2,2,2,2))
        }
      }
      par(mfrow=c(1,1))
      if(plot_options$save)
        dev.off()
      #  par(oldpar)
    }else{

      for(i in seq(length(CI_const))){
        if(plot_options$save)
          pdf(file = paste(plot_options$folderPlots,"profMeanCI",CI_const[i]*100,"_lower",plot_options$id_save,".pdf",sep=""),width = 18,height = 9)
        plotOneBivProfile(allRes = prof_CI_const_full[[i]]$lower,allPsi = allPsi,Design=object$profMean_full$Design,threshold=threshold,main_addendum=paste("(CI,",CI_const[i],", lower)"),xlab=expression(eta[1]),ylab=expression(eta[2]))
        if(plot_options$save)
          dev.off()

        if(plot_options$save)
          pdf(file = paste(plot_options$folderPlots,"profMeanCI",CI_const[i]*100,"_upper",plot_options$id_save,".pdf",sep=""),width = 18,height = 9)
        plotOneBivProfile(allRes = prof_CI_const_full[[i]]$upper,allPsi = allPsi,Design=object$profMean_full$Design,threshold=threshold,main_addendum=paste("(CI,",CI_const[i],", upper)"),xlab=expression(eta[1]),ylab=expression(eta[2]))
        if(plot_options$save)
          dev.off()
      }
    }
    object$profCI_full<-prof_CI_const_full
  }


  ## UQ part ##
  #  res_UQ<-NULL
  if(uq_computations){# && is.null(object$res_UQ)){
  #  warning("UQ computations need to be checked!")
    if(is.null(object$res_UQ$kmModel))
      object$res_UQ$kmModel=object$kmModel
    object$res_UQ<-obliqueProf_UQ(object=object$res_UQ,threshold=threshold,allPsi=allPsi,allResMean=object$profMean_full,
                                options_approx=options_approx,plot_level=max(0,plot_level-1),
                                plot_options=plot_options,return_level=max(1,return_level-1),...)

  }else{
    object$res_UQ<-NULL
  }
  #results<-object =list(kmModel=object$kmModel, profMean_full=object$profMean_full,profMean_approx=object$profMean_approx,res_UQ=object$res_UQ)
  if(return_level==1){
    return(object)
  }else{
    if(is.null(object$more)){
      times<-list(full=timeMean_full,approx=timeMean_approx)
      object$more<-list(abs_err=abs_err,times=times)
    }
    return(object)
  }

}
