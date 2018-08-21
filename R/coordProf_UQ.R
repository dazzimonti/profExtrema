# coordProf_UQ function
#' @name coordProf_UQ
#' @author Dario Azzimonti
#' @title Coordinate profiles UQ from a kriging model
#' @description The function coordProf_UQ computes the profile extrema functions for posterior realizations of a Gaussian process and its confidence bounds
#' @param object either a \link[DiceKriging]{km} model or a list containing partial results. If \code{object} is a km model then all computations are carried out. If \code{object} is a list, then the function carries out all computations to complete the results list.
#' @param threshold the threshold of interest
#' @param allResMean a list resulting from \code{getAllMaxMin} or \code{approxMaxMin} for the profile extrema on the mean. If NULL the median from the observations is plotted
#' @param quantiles_uq a vector containing the quantiles to be computed
#' @param options_approx an optional list of options for approxMaxMin, see \link{approxMaxMin} for details.
#' @param options_full_sims an optional list of options for getAllMaxMin, see \link{getAllMaxMin} for details. If NULL the full computations are not excuted. NOTE: this computations might be very expensive!
#' @param options_sims an optional list of options for the posterior simulations.
#' \itemize{
#' \item{\code{algorithm:}} string choice of the algorithm to select the pilot points ("A" or "B", default "B");
#' \item{\code{lower:}} \eqn{d} dimensional vector with lower bounds for pilot points, default \code{rep(0,d)};
#' \item{\code{upper:}} \eqn{d} dimensional vector with upper bounds for pilot points, default \code{rep(1,d)};
#' \item{\code{batchsize:}} number of pilot points, default \code{120};
#' \item{\code{optimcontrol:}} list containing the options for optimization, see \link[pGPx]{optim_dist_measure};
#' \item{\code{integcontrol:}} list containing the options for numerical integration of the criterion, see \link[pGPx]{optim_dist_measure};
#' \item{\code{integration.param:}} list containing the integration design, obtained with the function \link[KrigInv]{integration_design};
#' \item{\code{nsim:}} number of approximate GP simulations, default \code{300}.
#' }
#' @param options_bound an optional list containing \code{beta} the confidence level for the approximation and \code{alpha} the confidence level for the bound. Note that \code{alpha > 2*beta}. If \code{NULL}, the bound is not computed.
#' @param plot_level an integer to select the plots to return (0=no plots, 1=basic plots, 2= all plots)
#' @param plot_options an optional list of parameters for plots. See \link{setPlotOptions} for currently available options.
#' @param return_level an integer to select the amount of details returned
#' @return If return_level=1 a list containing \itemize{
#' \item{\code{profSups:}}{an array \code{dxfullDesignSizexnsims} containing the profile sup for each coordinate for each realization.}
#' \item{\code{profInfs:}}{an array \code{dxfullDesignSizexnsims} containing the profile inf for each coordinate for each realization.}
#' \item{\code{prof_quantiles_approx:}}{a list containing the quantiles (levels set by \code{quantiles_uq}) of the profile extrema functions.}
#' } if return_level=2 the same list as above but also including \code{more:} a list containing \itemize{
#' \item{\code{times:}}{a list containing
#' 	\itemize{
#' 	\item{\code{tSpts:} }{computational time for selecting pilot points.}
#' 	\item{\code{tApprox1ord:}}{vector containing the computational time required for profile extrema computation for each realization}
#' 	}}
#' \item{\code{simuls:}}{ a matrix containing the value of the field simulated at the pilot points}
#' \item{\code{sPts:}}{the pilot points}
#' }
#' @examples
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # Compute a kriging model from 50 evaluations of the Branin function
#' # Define the function
#' g<-function(x){
#'   return(-branin(x))
#' }
#' gp_des<-lhs::maximinLHS(20,2)
#' reals<-apply(gp_des,1,g)
#' kmModel<-km(design = gp_des,response = reals,covtype = "matern3_2")
#'
#' threshold=-10
#' d<-2
#'
#' # Compute coordinate profiles UQ starting from GP model
#' # define simulation options
#' options_sims<-list(algorithm="B", lower=rep(0,d), upper=rep(1,d),
#'                    batchsize=80, optimcontrol = list(method="genoud",pop.size=100,print.level=0),
#'                    integcontrol = list(distrib="sobol",n.points=1000), nsim=150)
#' # define 1 order approximation options
#' init_des<-lhs::maximinLHS(15,d)
#' options_approx<- list(multistart=4,heavyReturn=TRUE,
#'                       initDesign=init_des,fullDesignSize=100,
#'                       smoother="1order")
#' # define plot options
#' options_plots<-list(save=FALSE, titleProf = "Coordinate profiles",
#'                     title2d = "Posterior mean",qq_fill=TRUE)
#' \dontrun{
#' # profile UQ on approximate coordinate profiles
#' cProfiles_UQ<-coordProf_UQ(object = kmModel,threshold = threshold,allResMean = NULL,
#'                             quantiles_uq = c(0.05,0.95),options_approx = options_approx,
#'                             options_full_sims = NULL,options_sims = options_sims,
#'                             options_bound = NULL,plot_level = 3,
#'                             plot_options = options_plots,return_level = 3)
#' # profile UQ on full optim coordinate profiles
#' options_full_sims<-list(multistart=4,heavyReturn=TRUE)
#' cProfiles_UQ_full<-coordProf_UQ(object = cProfiles_UQ,threshold = threshold,allResMean = NULL,
#'                             quantiles_uq = c(0.05,0.95),options_approx = options_approx,
#'                             options_full_sims = options_full_sims,options_sims = options_sims,
#'                             options_bound = NULL,plot_level = 3,
#'                             plot_options = options_plots,return_level = 3)
#'
#' # profile UQ on full optim coordinate profiles with bound
#' cProfiles_UQ_full_bound<-coordProf_UQ(object = cProfiles_UQ_full,threshold = threshold,
#'                                       allResMean = NULL, quantiles_uq = c(0.05,0.95),
#'                                       options_approx = options_approx,
#'                                       options_full_sims = options_full_sims,
#'                                       options_sims = options_sims,
#'                                       options_bound = list(beta=0.024,alpha=0.05),
#'                                       plot_level = 3, plot_options = options_plots,
#'                                       return_level = 3)
#' }
#' @export
coordProf_UQ = function(object,threshold,allResMean=NULL,quantiles_uq=c(0.05,0.95),options_approx=NULL,options_full_sims=NULL,options_sims=NULL,options_bound=NULL,plot_level=0,plot_options=NULL,return_level=1){

  # number of thresholds
  num_T<-length(threshold)

  # Check object
  if(is(object,"km")){
    object<-list(kmModel=object)
  }else if(!is.list(object)){
    stop("object must be either a list or a km object")
  }

  # set up dimension
  d<-object$kmModel@d

  # Options setup
  if(is.null(options_approx)){
    init_des<-maximinLHS(10,2)
    options_approx<- list(multistart=8,heavyReturn=TRUE,initDesign=init_des,fullDesignSize=100)
  }

  # Set up plot options
  plot_options<-setPlotOptions(plot_options = plot_options,d=d,num_T=num_T,kmModel=object$kmModel)

  # Set-up simulation options
  if(is.null(options_sims)){
    options_sims<-list(algorithm="B", lower=rep(0,d), upper=rep(1,d),
                       batchsize=120, optimcontrol = list(method="genoud",pop.size=100,print.level=0),
                       integcontrol = list(distrib="sobol",n.points=1000),
                       nsim=300)
    options_sims$integration.param = integration_design(options_sims$integcontrol,d,options_sims$lower,options_sims$upper,object$kmModel,threshold)
    options_sims$integration.param$alpha <- 0.5
  }

  if(is.null(allResMean)){
    quantiles_uq<-c(quantiles_uq,0.5)
  }else{
    changePP<-getChangePoints(threshold = threshold,allRes = allResMean)
  }
  ##### Get the pilot points


  # If not already in object, obtain the pilot points
  if(is.null(object$sPts)){
    timeIn<-get_nanotime()
    object$sPts<-optim_dist_measure(model = object$kmModel,threshold = threshold[1],
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

  ###
  # Prepare the functions for UQ profiles
  nugget.sim=1e-6
  type="UK"
  simu_points<-object$sPts$par

  if(is.null(object$more$simuls)){
    some.simu <- simulate(object=object$kmModel,nsim=options_sims$nsim,newdata=simu_points,nugget.sim=nugget.sim,
                             cond=TRUE,checkNames = FALSE)
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
    object$profSups<-array(NA,dim = c(d,options_approx$fullDesignSize,options_sims$nsim))
    object$profInfs<-array(NA,dim = c(d,options_approx$fullDesignSize,options_sims$nsim))
    tApprox1ord<-rep(NA,options_sims$nsim)
  }

  if(!is.null(options_full_sims) && is.null(object$profSups_full)){
    object$profSups_full<-array(NA,dim = c(d,options_approx$fullDesignSize,options_sims$nsim))
    object$profInfs_full<-array(NA,dim = c(d,options_approx$fullDesignSize,options_sims$nsim))
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
      temp_full<-getAllMaxMin(f = g_uq_spec,fprime = g_uq_der_spec,d = d,options = options_full_sims)
      tFull[i]<-(get_nanotime()-timeIn)*1e-9

      object$profSups_full[,,i]<-t(temp_full$res$max)
      object$profInfs_full[,,i]<-t(temp_full$res$min)
    }

    if(is.na(object$profSups[1,1,i]) || is.na(object$profInfs[1,1,i])){
      if(i%%10==0){
        cat("Approx_sims. Realization ",i,"\n")
      }
      timeIn<-get_nanotime()
      temp_1o<-approxMaxMin(f = g_uq_spec,fprime = g_uq_der_spec,d = d,opts = options_approx)
      tApprox1ord[i]<-(get_nanotime()-timeIn)*1e-9

      #  temp<-getAllMaxMin(f=g_uq_spec,fprime = NULL,d=2,options = list(multistart=2,heavyReturn=TRUE))
      object$profSups[,,i]<-t(temp_1o$res$max)
      object$profInfs[,,i]<-t(temp_1o$res$min)
    }
  }

  # save quantiles for approximations
  object$prof_quantiles_approx<-list()
  for(i in seq(length(quantiles_uq))){
    object$prof_quantiles_approx[[i]]<-list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d),
                                                     max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d)))
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
      object$prof_quantiles_full[[i]]<-list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d),
                                                     max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d)))
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
    image(matrix(pred2d$mean,nrow = 100),col=gray.colors(20), main=plot_options$title2d,xlab = plot_options$coord_names[1],ylab= plot_options$coord_names[2],
          cex.main=3,cex.axis=1.8,cex.lab=2.8)
    contour(matrix(pred2d$mean,nrow = 100),add=T,nlevels = 10,lwd=1.5,labcex=1.2)
    contour(matrix(pred2d$mean,nrow = 100),add=T,levels = threshold,col=plot_options$col_thresh,lwd=3,labcex=1.5)
    for(tt in seq(num_T)){
      abline(v = changePP$neverEx[[tt]][[1]],col=plot_options$col_CCPthresh_nev[tt],lwd=2.5)
      abline(h = changePP$neverEx[[tt]][[2]],col=plot_options$col_CCPthresh_nev[tt],lwd=2.5)
      abline(v = changePP$alwaysEx[[tt]][[1]],col=plot_options$col_CCPthresh_alw[tt],lwd=2.5)
      abline(h = changePP$alwaysEx[[tt]][[2]],col=plot_options$col_CCPthresh_alw[tt],lwd=2.5)
    }
    if(!is.null(options_full_sims)){
      for(j in seq(length(quantiles_uq))){
        for(tt in seq(num_T)){
          abline(v = ccPP_full[[j]]$neverEx[[tt]][[1]],col=plot_options$col_CCPthresh_nev[tt],lwd=2,lty=2)
          abline(h = ccPP_full[[j]]$neverEx[[tt]][[2]],col=plot_options$col_CCPthresh_nev[tt],lwd=2,lty=2)
          abline(v = ccPP_full[[j]]$alwaysEx[[tt]][[1]],col=plot_options$col_CCPthresh_alw[tt],lwd=2,lty=2)
          abline(h = ccPP_full[[j]]$alwaysEx[[tt]][[2]],col=plot_options$col_CCPthresh_alw[tt],lwd=2,lty=2)
        }
      }
    }else{
      for(j in seq(length(quantiles_uq))){
        for(tt in seq(num_T)){
          abline(v = ccPP[[j]]$neverEx[[tt]][[1]],col=plot_options$col_CCPthresh_nev[tt],lwd=2,lty=2)
          abline(h = ccPP[[j]]$neverEx[[tt]][[2]],col=plot_options$col_CCPthresh_nev[tt],lwd=2,lty=2)
          abline(v = ccPP[[j]]$alwaysEx[[tt]][[1]],col=plot_options$col_CCPthresh_alw[tt],lwd=2,lty=2)
          abline(h = ccPP[[j]]$alwaysEx[[tt]][[2]],col=plot_options$col_CCPthresh_alw[tt],lwd=2,lty=2)
        }
      }
    }
    if(plot_options$fun_evals>0){
      points(object$kmModel@X,pch=17,cex=1.6)
    }
    if(plot_options$save)
      dev.off()
  }

  if(plot_level>=1){
    object$bound$bound <- NULL
    plot_univariate_profiles_UQ(objectUQ = object, plot_options = plot_options,nsims = options_sims$nsim,quantiles_uq=quantiles_uq,
                                threshold = threshold,nameFile ="prof_UQ_approx", profMean = allResMean,typeProf = "approx")

    if(!is.null(options_full_sims))
      plot_univariate_profiles_UQ(objectUQ = object, plot_options = plot_options,nsims = options_sims$nsim,quantiles_uq=quantiles_uq,
                                  threshold = threshold,nameFile ="prof_UQ_full", profMean = allResMean,typeProf = "full")
  }
  #  object$profSups=profSups
  #  object$profInfs=profInfs
  #  object$prof_quantiles_approx=prof_quantiles_approx
  #  object$sPts=m_dist

  # Compute the bound correction
  if(!is.null(options_bound)){
    object$bound<-bound_profiles(objectUQ = object,mean_var_delta = object$bound$mean_var_D,beta = options_bound$beta,alpha = options_bound$alpha,
                                 options_approx = options_approx,options_full_sims = options_full_sims)

    if(plot_level>=1){

      plot_univariate_profiles_UQ(objectUQ = object, plot_options = plot_options,nsims = options_sims$nsim,quantiles_uq=quantiles_uq,
                                  threshold = threshold,nameFile ="prof_UQ_bound_approx", profMean = allResMean,typeProf = "approx")

      if(!is.null(options_full_sims))
        plot_univariate_profiles_UQ(objectUQ = object, plot_options = plot_options,nsims = options_sims$nsim,quantiles_uq=quantiles_uq,
                                    threshold = threshold,nameFile ="prof_UQ_bound_full", profMean = allResMean,typeProf = "full")
    }
  }


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
