# coordProf_UQ function
#' @author Dario Azzimonti
#' @title Coordinate profiles UQ from a kriging model
#' @description The function coordProf_UQ computes the profile extrema functions for posterior realizations of a Gaussian process and its confidence bounds
#' @param object either a \link[DiceKriging]{km} model or a list containing partial results. If \code{object} is a km model then all computations are carried out. If \code{object} is a list, then the function carries out all computations to complete the results list.
# @param kmModel a \link[DiceKriging]{km} model
#' @param threshold the threshold of interest
#' @param allResMean a list resulting from \code{getAllMaxMin} or \code{approxMaxMin} for the profile extrema on the mean. If NULL the median from the observations is plotted
#' @param quantiles_uq a vector containing the quantiles to be computed
#' @param options_approx an optional list of options for approxMaxMin, see \link{approxMaxMin} for details.
#' @param options_sims an optional list of options for the posterior simulations.
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
#' @return If return_level=1 a list containing \itemize{
#' \item{\code{profSups:}}{an array \code{dxfullDesignSizexnsims} containing the profile sup for each coordinate for each realization.}
#' \item{\code{profInfs:}}{an array \code{dxfullDesignSizexnsims} containing the profile inf for each coordinate for each realization.}
#' \item{\code{prof_quantiles_approx:}}{a list containing the quantiles (levels set by \code{quantiles_uq}) of the profile extrema functions.}
#' } if return_level=2 the same list as above but also including \code{more:} a list containing \itemize{
#' \item{\code{times:}}{a list containing
#' 	\itemize{
#' 	\item{\code{tSpts:} }{computational time for selecting simulation points.}
#' 	\item{\code{tApprox1ord:}}{vector containing the computational time required for profile extrema computation for each realization}
#' 	}}
#' \item{\code{simuls:}}{ a matrix containing the value of the field simulated at the simulation points}
#' \item{\code{sPts:}}{the simulation points}
#' }
#' @export
coordProf_UQ = function(object,threshold,allResMean=NULL,quantiles_uq=c(0.05,0.95),options_approx=NULL,options_sims=NULL,plot_level=0,plot_options=NULL,return_level=1){

  # Check object
  if(is(object,"km")){
    object<-list(kmModel=object)
  }else if(!is.list(object)){
    stop("object must be either a list or a km object")
  }

  # Options setup
  if(is.null(options_approx)){
    init_des<-maximinLHS(10,2)
    options_approx<- list(multistart=8,heavyReturn=TRUE,initDesign=init_des,fullDesignSize=100)
  }

  if(is.null(plot_options)){
    if(plot_level>0){
      plot_options<-list(save=F)
    }

  }else{
    if(plot_options$save==TRUE && is.null(plot_options$folderPlots))
      plot_options$folderPlots <- './'

    if(is.null(plot_options$design)){
      plot_options$design<-matrix(NA,ncol=object$kmModel@d,nrow=100)
      for(i in seq(object$kmModel@d)){
        plot_options$design[,i]<-seq(0,1,,100)
      }
    }
    # Useful for serial plots in HPC
    if(is.null(plot_options$id_save)){
      plot_options$id_save<-""
    }

    # if true it fills the region between the first 2 quantiles in quantiles_uq
    if(is.null(plot_options$qq_fill)){
      plot_options$qq_fill<-FALSE
    }
  }

  if(is.null(options_sims)){
    options_sims<-list(algorithm="B", lower=rep(0,object$kmModel@d), upper=rep(1,object$kmModel@d),
                       batchsize=120, optimcontrol = list(method="genoud",pop.size=100,print.level=0),
                       integcontrol = list(distrib="sobol",n.points=1000),
                       nsim=300)
    options_sims$integration.param = integration_design(options_sims$integcontrol,object$kmModel@d,options_sims$lower,options_sims$upper,object$kmModel,threshold)
    options_sims$integration.param$alpha <- 0.5
  }

  if(is.null(allResMean)){
    quantiles_uq<-c(quantiles_uq,0.5)
  }else{
    changePP<-getChangePoints(threshold = threshold,allRes = allResMean)
  }
  ##### Get the simulation points


  # If not already in object, obtain the simulation points
  if(is.null(object$sPts)){
    timeIn<-get_nanotime()
    object$sPts<-optim_dist_measure(model = object$kmModel,threshold = threshold,
                                    lower = options_sims$lower,upper = options_sims$upper,batchsize = options_sims$batchsize,
                                    algorithm = options_sims$algorithm,verb=1,optimcontrol = options_sims$optimcontrol,integration.param = options_sims$integration.param)
    timeMdist<-(get_nanotime()-timeIn)*1e-9
  }else{
    if(is.null(object$more$times$tSpts)){
      times<-list(tSpts=NA)
      if(is.null(object$more)){
        object$more<-list(times=times)
      }else{
        object$more$times=times
      }
    }
  }


  #cairo_pdf(paste(plot_options$folderPlots,"critVal_tt",index_exp,".pdf",sep=""),width = 14,height = 14)
  if(plot_level>0){
    oldpar<-par()
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"sPtsCritVal",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    plot(object$sPts$value,type='o',main="Optimized criterion",ylab="y",xlab="Iteration")
    if(plot_options$save)
      dev.off()
  }

  d<-object$kmModel@d

  ###
  # if the profSups and profInfs are not already there, compute them
  if(is.null(object$profSups) || is.null(object$profInfs)){
    nugget.sim=1e-6
    type="UK"
    simu_points<-object$sPts$par

    if(is.null(object$more$simuls)){
      some.simu <- simulate_km(object=object$kmModel,nsim=options_sims$nsim,newdata=simu_points,nugget.sim=nugget.sim,
                               cond=TRUE,checkNames = FALSE, type=type)
    }else{
      some.simu<-object$more$simuls
    }

    g_uq<-function(x,realization,kmModel,simupoints,F.mat=NULL,T.mat=NULL){
      x<-matrix(x,ncol=kmModel@d)
      colnames(x)<-colnames(kmModel@X)
      obj <- krig_weight_GPsimu(object=kmModel,simu_points=simupoints,krig_points=x,T.mat = T.mat,F.mat = F.mat)
      krig.mean.init <- matrix(obj$krig.mean.init,ncol=1)
      weights <- t(obj$Lambda.end)

      return(krig.mean.init + tcrossprod(weights,matrix(realization,nrow=1)))
    }

    g_uq_deriv<-function(x,realization,kmModel,simupoints,T.mat=NULL,F.mat=NULL){
      x<-matrix(x,ncol=kmModel@d)
      colnames(x)<-colnames(kmModel@X)
      obj_deriv<-grad_kweights(object = kmModel,simu_points = simupoints,krig_points = matrix(x,ncol=kmModel@d),T.mat = T.mat,F.mat = F.mat)
      krig_mean_init <- matrix(obj_deriv$krig.mean.init,ncol=kmModel@d)
      weights <- t(obj_deriv$Lambda.end)

      return(krig_mean_init + tcrossprod(matrix(realization,nrow=1),weights))
    }


    ### Lets compute the profile extrema for this realization
    # choose size of full design
    object$profSups<-array(NA,dim = c(object$kmModel@d,options_approx$fullDesignSize,options_sims$nsim))
    object$profInfs<-array(NA,dim = c(object$kmModel@d,options_approx$fullDesignSize,options_sims$nsim))
    tApprox1ord<-rep(NA,options_sims$nsim)


    F.mat <- model.matrix(object=object$kmModel@trend.formula, data = data.frame(rbind(object$kmModel@X,simu_points)))
    K <- covMatrix(object=object$kmModel@covariance,X=rbind(object$kmModel@X,simu_points))$C
    T.mat <- chol(K)
    for(i in seq(options_sims$nsim)){

      g_uq_spec<-function(x){
        return(g_uq(x=x,realization=some.simu[i,],kmModel = object$kmModel,simupoints = simu_points,F.mat = F.mat,T.mat = T.mat))
      }
      g_uq_der_spec<-function(x){
        return(g_uq_deriv(x=x,realization=some.simu[i,],kmModel = object$kmModel,simupoints = simu_points,T.mat = T.mat,F.mat = F.mat))
      }
      if(i%%10==0){
        cat("Realization ",i,"\n")
      }

      timeIn<-get_nanotime()
      temp_1o<-approxMaxMin(f = g_uq_spec,fprime = g_uq_der_spec,threshold = threshold,d = object$kmModel@d,opts = options_approx)
      tApprox1ord[i]<-(get_nanotime()-timeIn)*1e-9

      #  temp<-getAllMaxMin(f=g_uq_spec,fprime = NULL,d=2,options = list(multistart=2,heavyReturn=TRUE))
      object$profSups[,,i]<-t(temp_1o$res$max)
      object$profInfs[,,i]<-t(temp_1o$res$min)
    }
  }


  object$prof_quantiles_approx<-list()
  for(i in seq(length(quantiles_uq))){
    object$prof_quantiles_approx[[i]]<-list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = object$kmModel@d),
                                                     max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = object$kmModel@d)))
  }
  names(object$prof_quantiles_approx)<-quantiles_uq

  #qq95_sup<-matrix(NA,nrow = 100,ncol = 2)
  #qq05_sup<-matrix(NA,nrow = 100,ncol = 2)
  #qq95_inf<-matrix(NA,nrow = 100,ncol = 2)
  #qq05_inf<-matrix(NA,nrow = 100,ncol = 2)

  ccPP<-list()
  for(j in seq(length(quantiles_uq))){
    for(coord in seq(d)){
      object$prof_quantiles_approx[[j]]$res$max[,coord]<-apply(object$profSups[coord,,],1,function(x){return(quantile(x,quantiles_uq[j]))})
      object$prof_quantiles_approx[[j]]$res$min[,coord]<-apply(object$profInfs[coord,,],1,function(x){return(quantile(x,quantiles_uq[j]))})
    }
    ccPP[[j]]<-getChangePoints(threshold = threshold,allRes = object$prof_quantiles_approx[[j]])
  }
  names(ccPP)<-quantiles_uq

  ## Plot profiles with Uncertainty
  dd<-seq(0,1,,length.out = options_approx$fullDesignSize)

  if(is.null(allResMean))
    changePP<-ccPP$`0.5`

  # Plot the posterior mean and visualize the actual excursion set and the regions of no-excursion according to the profile extrema functions.
  if(plot_level>=2 && d==2){
    # since dimension==2 we can plot the posterior mean
    newdata<-expand.grid(seq(0,1,,100),seq(0,1,,100))
    colnames(newdata)<-colnames(object$kmModel@X)
    pred2d<-predict.km(object$kmModel,newdata = newdata,type = "UK",light.return = TRUE,se.compute = FALSE)

    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"profMean_UQ",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    par(mar = c(5, 5, 4, 2) + 0.1)
    image(matrix(pred2d$mean,nrow = 100),col=gray.colors(20), main=plot_options$title2d,xlab = colnames(object$kmModel@X)[1],ylab= colnames(object$kmModel@X)[2],
          cex.main=3,cex.axis=1.8,cex.lab=2.8)
    contour(matrix(pred2d$mean,nrow = 100),add=T,nlevels = 10,lwd=1.5,labcex=1.2)
    contour(matrix(pred2d$mean,nrow = 100),add=T,levels = threshold,col=2,lwd=3,labcex=1.5)
    abline(v = changePP$neverEx[[1]],col=3,lwd=2.5)
    abline(h = changePP$neverEx[[2]],col=3,lwd=2.5)
    for(j in seq(length(quantiles_uq))){
      abline(v = ccPP[[j]]$neverEx[[1]],col=4,lwd=2,lty=2)
      abline(h = ccPP[[j]]$neverEx[[2]],col=4,lwd=2,lty=2)
      abline(v = ccPP[[j]]$alwaysEx[[1]],col=5,lwd=2,lty=2)
      abline(h = ccPP[[j]]$alwaysEx[[2]],col=5,lwd=2,lty=2)
    }
    if(plot_options$save)
      dev.off()
  }

  if(plot_level>=1){
    if(plot_options$save)
      cairo_pdf(filename = paste(plot_options$folderPlots,"prof_UQ",plot_options$id_save,".pdf",sep=""),width = 12,height = 12)
    mfrows<-switch(d,"1"=c(1,1),"2"=c(2,1),"3"=c(2,2),"4"=c(2,2),"5"=c(3,2),"6"=c(3,2))
    par(mfrow=mfrows, mar = c(4, 5, 3, 1) + 0.1)
    for(coord in seq(d)){
      title_string<-paste("Coordinate",colnames(object$kmModel@X)[coord])
      if(!is.null(allResMean)){
        ylimTemp<-range(c(range(object$profSups[coord,,]),range(object$profInfs[coord,,]),
                          range(allResMean$res$min[,coord]),range(allResMean$res$max[,coord])))
        #ylimTemp<-c(min(allResMean$res$min[,coord]),max(allResMean$res$max[,coord]))
        plot(plot_options$design[,coord],allResMean$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      }else{
        ylimTemp<-range(c(range(object$profSups[coord,,]),range(object$profInfs[coord,,]),
                          range(object$prof_quantiles_approx$`0.5`$res$min[,coord]),
                          range(object$prof_quantiles_approx$`0.5`$res$max[,coord])))
        plot(plot_options$design[,coord],object$prof_quantiles_approx$`0.5`$res$min[,coord],ylim=ylimTemp,type='l',main=title_string,
             xlab=expression(eta),ylab="f",lty=1,cex.main=3,cex.lab=2.5,cex.axis=1.8)
      }
      for(i in seq(options_sims$nsim)){
        dd<-seq(0,1,,length.out = length(object$profInfs[coord,,i]))
        lines(dd, object$profSups[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
        lines(dd, object$profInfs[coord,,i],col=adjustcolor('darkolivegreen3',alpha.f=0.5),lwd=0.2)
      }
      if(plot_options$qq_fill){
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(object$prof_quantiles_approx[[2]]$res$min[,coord],
                  rev(object$prof_quantiles_approx[[1]]$res$min[,coord])),
                col = adjustcolor("red",alpha.f=0.5), border = NA)
        polygon(c(plot_options$design[,coord], rev(plot_options$design[,coord])),
                c(object$prof_quantiles_approx[[2]]$res$max[,coord],
                  rev(object$prof_quantiles_approx[[1]]$res$max[,coord])),
                col = adjustcolor("red",alpha.f=0.5), border = NA,alpha=0.5)
      }
      #  lines(plot_options$design[,coord],allRes_approx$res$max[,coord],lty=3,col=4,lwd=2)
      #  points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$max[,coord],col=4)
      #  lines(plot_options$design[,coord],allRes_approx$res$min[,coord],lty=4,col=4,lwd=2)
      #  points(allRes_approx$profPoints$design[,coord],allRes_approx$profPoints$res$min[,coord],col=4)
      abline(v=changePP$neverEx[[coord]],col=3,lwd=2.5)
      abline(h = threshold,col=2,lwd=2)
      for(j in seq(length(quantiles_uq))){
        lines(plot_options$design[,coord],object$prof_quantiles_approx[[j]]$res$max[,coord],lty=j,col=1,lwd=1.5)
        lines(plot_options$design[,coord],object$prof_quantiles_approx[[j]]$res$min[,coord],lty=j,col=1,lwd=1.5)
        abline(v=ccPP[[j]]$neverEx[[coord]],col=4,lwd=2,lty=2)
        abline(v=ccPP[[j]]$neverEx[[coord]],col=4,lwd=2,lty=3)
      }
      if(!is.null(allResMean)){
        lines(plot_options$design[,coord],allResMean$res$max[,coord],lty=1)
        lines(plot_options$design[,coord],allResMean$res$min[,coord],lty=1)
      }else{
        lines(plot_options$design[,coord],object$prof_quantiles_approx$`0.5`$res$max[,coord],lty=1)
        lines(plot_options$design[,coord],object$prof_quantiles_approx$`0.5`$res$min[,coord],lty=1)
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
    #  par(oldpar)
  }
  #  object$profSups=profSups
  #  object$profInfs=profInfs
  #  object$prof_quantiles_approx=prof_quantiles_approx
  #  object$sPts=m_dist

  if(return_level==1){
    return(object)
  }else{
    if(is.null(object$more)){
      times<-list(tSpts=timeMdist,tApprox1ord=tApprox1ord)
      object$more<-list(simuls=some.simu,times=times)
    }
    return(object)
  }

}
