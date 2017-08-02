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
#' @param options_full_sims an optional list of options for getAllMaxMin, see \link{getAllMaxMin} for details. If NULL the full computations are not excuted. NOTE: this computations might be very expensive!
#' @param options_sims an optional list of options for the posterior simulations.
#' \itemize{
#' \item{\code{algorithm:}} string choice of the algorithm to select the simulation points ("A" or "B");
#' \item{\code{lower:}} \eqn{d} dimensional vector with lower bounds for simulation points;
#' \item{\code{upper:}} \eqn{d} dimensional vector with upper bounds for simulation points;
#' \item{\code{batchsize:}} number of simulation points;
#' \item{\code{optimcontrol:}} list containing the options for optimization;
#' \item{\code{integcontrol:}} list containing the options for numerical integration of the criterion;
#' \item{\code{integration.param:}} list containing the integration design, obtained with the function \link[KrigInv]{integration_design}.
#' }
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
coordProf_UQ = function(object,threshold,allResMean=NULL,quantiles_uq=c(0.05,0.95),options_approx=NULL,options_full_sims=NULL,options_sims=NULL,plot_level=0,plot_options=NULL,return_level=1){

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
  # Prepare the functions for UQ profiles
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

  # Useful one-time computations
  F.mat <- model.matrix(object=object$kmModel@trend.formula, data = data.frame(rbind(object$kmModel@X,simu_points)))
  K <- covMatrix(object=object$kmModel@covariance,X=rbind(object$kmModel@X,simu_points))$C
  T.mat <- chol(K)

  ### Lets compute the profile extrema for this realization
  # if the profSups and profInfs are not already there, compute them
  if(is.null(object$profSups) || is.null(object$profInfs)){
    # choose size of full design
    object$profSups<-array(NA,dim = c(object$kmModel@d,options_approx$fullDesignSize,options_sims$nsim))
    object$profInfs<-array(NA,dim = c(object$kmModel@d,options_approx$fullDesignSize,options_sims$nsim))
    tApprox1ord<-rep(NA,options_sims$nsim)
  }

  if(!is.null(options_full_sims) && is.null(object$profSups_full)){
    object$profSups_full<-array(NA,dim = c(object$kmModel@d,options_approx$fullDesignSize,options_sims$nsim))
    object$profInfs_full<-array(NA,dim = c(object$kmModel@d,options_approx$fullDesignSize,options_sims$nsim))
    tFull<-rep(NA,options_sims$nsim)
  }

  for(i in seq(options_sims$nsim)){

    g_uq_spec<-function(x){
      return(g_uq(x=x,realization=some.simu[i,],kmModel = object$kmModel,simupoints = simu_points,F.mat = F.mat,T.mat = T.mat))
    }
    g_uq_der_spec<-function(x){
      return(g_uq_deriv(x=x,realization=some.simu[i,],kmModel = object$kmModel,simupoints = simu_points,T.mat = T.mat,F.mat = F.mat))
    }


    if(!is.null(options_full_sims) && is.na(object$profSups_full[1,1,i])){
      if(i%%10==0){
        cat("Full_sims.Realization ",i,"\n")
      }
      timeIn<-get_nanotime()
      temp_full<-getAllMaxMin(f = g_uq_spec,fprime = g_uq_der_spec,d = object$kmModel@d,options = options_full_sims)
      tFull[i]<-(get_nanotime()-timeIn)*1e-9

      object$profSups_full[,,i]<-t(temp_full$res$max)
      object$profInfs_full[,,i]<-t(temp_full$res$min)
    }

    if(is.null(object$profSups) || is.null(object$profInfs)){
      if(i%%10==0){
        cat("Approx_sims. Realization ",i,"\n")
      }
      timeIn<-get_nanotime()
      temp_1o<-approxMaxMin(f = g_uq_spec,fprime = g_uq_der_spec,threshold = threshold,d = object$kmModel@d,opts = options_approx)
      tApprox1ord[i]<-(get_nanotime()-timeIn)*1e-9

      #  temp<-getAllMaxMin(f=g_uq_spec,fprime = NULL,d=2,options = list(multistart=2,heavyReturn=TRUE))
      object$profSups[,,i]<-t(temp_1o$res$max)
      object$profInfs[,,i]<-t(temp_1o$res$min)
    }
  }

  # save quantiles for approximations
  object$prof_quantiles_approx<-list()
  for(i in seq(length(quantiles_uq))){
    object$prof_quantiles_approx[[i]]<-list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = object$kmModel@d),
                                                     max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = object$kmModel@d)))
  }
  names(object$prof_quantiles_approx)<-quantiles_uq

  ccPP<-list()
  for(j in seq(length(quantiles_uq))){
    for(coord in seq(d)){
      object$prof_quantiles_approx[[j]]$res$max[,coord]<-apply(object$profSups[coord,,],1,function(x){return(quantile(x,quantiles_uq[j]))})
      object$prof_quantiles_approx[[j]]$res$min[,coord]<-apply(object$profInfs[coord,,],1,function(x){return(quantile(x,quantiles_uq[j]))})
    }
    ccPP[[j]]<-getChangePoints(threshold = threshold,allRes = object$prof_quantiles_approx[[j]])
  }
  names(ccPP)<-quantiles_uq

  # save quantiles for full optim
  if(!is.null(options_full_sims)){
    object$prof_quantiles_full<-list()
    for(i in seq(length(quantiles_uq))){
      object$prof_quantiles_full[[i]]<-list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = object$kmModel@d),
                                                     max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = object$kmModel@d)))
    }
    names(object$prof_quantiles_full)<-quantiles_uq

    ccPP_full<-list()
    for(j in seq(length(quantiles_uq))){
      for(coord in seq(d)){
        object$prof_quantiles_full[[j]]$res$max[,coord]<-apply(object$profSups_full[coord,,],1,function(x){return(quantile(x,quantiles_uq[j]))})
        object$prof_quantiles_full[[j]]$res$min[,coord]<-apply(object$profInfs_full[coord,,],1,function(x){return(quantile(x,quantiles_uq[j]))})
      }
      ccPP_full[[j]]<-getChangePoints(threshold = threshold,allRes = object$prof_quantiles_full[[j]])
    }
    names(ccPP_full)<-quantiles_uq

  }

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

    plot_univariate_profiles_UQ(objectUQ = object, plot_options = plot_options,nsims = options_sims$nsim,
                                threshold = threshold,nameFile ="prof_UQ_approx", profMean = allResMean,typeProf = "approx")

    if(!is.null(options_full_sims))
      plot_univariate_profiles_UQ(objectUQ = object, plot_options = plot_options,nsims = options_sims$nsim,
                                threshold = threshold,nameFile ="prof_UQ_full", profMean = allResMean,typeProf = "full")
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
      if(!is.null(options_full_sims)){
        times$tFull<-tFull
      }
      object$more<-list(simuls=some.simu,times=times)
    }
    return(object)
  }

}
