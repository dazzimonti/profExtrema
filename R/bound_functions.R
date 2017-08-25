# mean_Delta_T function
#' @author Dario Azzimonti
#' @name mean_Delta_T
#' @title mean function of difference process
#' @description The function mean_Delta_T computes the mean function of the difference process \eqn{Z_x - \widetilde{Z}_x} at \code{x}.
#' @param x a matrix \eqn{r x d} containing the \eqn{r} points where the function is to be computed.
#' @param kmModel the \link[DiceKriging]{km} model of the Gaussian process \eqn{Z}.
#' @param simupoints the matrix \eqn{l x d} containing the simulation points \eqn{E}.
#' @param T.mat the upper triangular factor of the Choleski decomposition of the covariance matrix of \code{rbind(kmModel@X,simupoints)}
#' @param F.mat the evaluation of the trend function at \code{rbind(kmModel@X,simupoints)}, see \link[stats]{model.matrix}.
#' @return the value of the mean function at \code{x} for the difference process \eqn{Z^\Delta = Z_x - \widetilde{Z}_x}.
mean_Delta_T<-function(x,kmModel,simupoints,T.mat,F.mat){
  x<-matrix(x,ncol=kmModel@d)
  colnames(x)<-colnames(kmModel@X)
  obj <- krig_weight_GPsimu(object=kmModel,simu_points=simupoints,krig_points=x,T.mat = T.mat,F.mat = F.mat)
  mu_n<-predict.km(object = kmModel,newdata = x,type = "UK",se.compute = FALSE,light.return = TRUE)$mean
  mu_E<-predict.km(object = kmModel,newdata = simupoints,type = "UK",se.compute = FALSE,light.return = TRUE)$mean
  return(mu_n-obj$krig.mean.init-crossprod(obj$Lambda.end,mu_E))
#  return(-tcrossprod(weights,predict.km(object = kmModel,newdata = simupoints,type = "UK",se.compute = FALSE,cov.compute = FALSE,light.return = FALSE)$mean))
}

# grad_mean_Delta_T function
#' @author Dario Azzimonti
#' @name grad_mean_Delta_T
#' @title Gradient of the mean function of difference process
#' @description The function grad_mean_Delta_T computes the gradient for the mean function of the difference process \eqn{Z_x - \widetilde{Z}_x} at \code{x}.
#' @param x a matrix \eqn{r x d} containing the \eqn{r} points where the function is to be computed.
#' @param kmModel the \link[DiceKriging]{km} model of the Gaussian process \eqn{Z}.
#' @param simupoints the matrix \eqn{l x d} containing the simulation points \eqn{E}.
#' @param T.mat the upper triangular factor of the Choleski decomposition of the covariance matrix of \code{rbind(kmModel@X,simupoints)}
#' @param F.mat the evaluation of the trend function at \code{rbind(kmModel@X,simupoints)}, see \link[stats]{model.matrix}.
#' @return the value of the gradient for the mean function at \code{x} for the difference process \eqn{Z^\Delta = Z_x - \widetilde{Z}_x}.
grad_mean_Delta_T<-function(x,kmModel,simupoints,T.mat,F.mat){
  x<-matrix(x,ncol=kmModel@d)
  colnames(x)<-colnames(kmModel@X)
  obj_deriv<-grad_kweights(object = kmModel,simu_points = simupoints,krig_points = matrix(x,ncol=kmModel@d),T.mat = T.mat,F.mat = F.mat)
  krig_mean_init <- matrix(obj_deriv$krig.mean.init,ncol=kmModel@d)
  weights <- t(obj_deriv$Lambda.end)

  krig_mean_init <- matrix(obj_deriv$krig.mean.init,ncol=kmModel@d)
  weights <- t(obj_deriv$Lambda.end)

  rr<-gradKm_dnewdata(object = kmModel,newdata=x,type = "UK",se.compute = F)$mean-krig_mean_init - tcrossprod(predict.km(object = kmModel,newdata = simupoints,type = "UK",se.compute = FALSE,cov.compute = FALSE,light.return = FALSE)$mean,weights)
  return(rr) # return(- tcrossprod(predict.km(object = kmModel,newdata = simupoints,type = "UK",se.compute = FALSE,cov.compute = FALSE,light.return = FALSE)$mean,weights))
}


# var_Delta_T function
#' @author Dario Azzimonti
#' @name var_Delta_T
#' @title Variance function of difference process
#' @description The function var_Delta_T computes the gradient for the variance function of the difference process \eqn{Z_x - \widetilde{Z}_x} at \code{x}.
#' @param x a matrix \eqn{r x d} containing the \eqn{r} points where the function is to be computed.
#' @param kmModel the \link[DiceKriging]{km} model of the Gaussian process \eqn{Z}.
#' @param simupoints the matrix \eqn{l x d} containing the simulation points \eqn{E}.
#' @param T.mat the upper triangular factor of the Choleski decomposition of the covariance matrix of \code{rbind(kmModel@X,simupoints)}
#' @param F.mat the evaluation of the trend function at \code{rbind(kmModel@X,simupoints)}, see \link[stats]{model.matrix}.
#' @return the value of the variance function at \code{x} for the difference process \eqn{Z^\Delta = Z_x - \widetilde{Z}_x}.
var_Delta_T<-function(x,kmModel,simupoints,T.mat,F.mat){
  x<-matrix(x,ncol=kmModel@d)
  colnames(x)<-colnames(kmModel@X)
  obj <- krig_weight_GPsimu(object=kmModel,simu_points=simupoints,krig_points=x,T.mat = T.mat,F.mat = F.mat)
#  krig.mean.init <- matrix(obj$krig.mean.init,nrow=1)
  weights <- t(obj$Lambda.end)
  var_post<-predict.km(object = kmModel,newdata = x,type = "UK",se.compute =TRUE)$sd^2
  kEE<-predict.km(object = kmModel,newdata = simupoints,type = "UK",cov.compute = TRUE,light.return = TRUE)$cov #covMatrix(object = kmModel@covariance,X = rbind(kmModel@X,simupoints))$C
  return(list(vv=var_post, vv_E= -weights%*%tcrossprod(kEE,weights)))
}


# grad_var_Delta_T function
#' @author Dario Azzimonti
#' @name grad_var_Delta_T
#' @title Gradient of the variance function of difference process
#' @description The function grad_var_Delta_T computes the gradient for the variance function of the difference process \eqn{Z_x - \widetilde{Z}_x} at \code{x}.
#' @param x a matrix \eqn{r x d} containing the \eqn{r} points where the function is to be computed.
#' @param kmModel the \link[DiceKriging]{km} model of the Gaussian process \eqn{Z}.
#' @param simupoints the matrix \eqn{l x d} containing the simulation points \eqn{E}.
#' @param T.mat the upper triangular factor of the Choleski decomposition of the covariance matrix of \code{rbind(kmModel@X,simupoints)}
#' @param F.mat the evaluation of the trend function at \code{rbind(kmModel@X,simupoints)}, see \link[stats]{model.matrix}.
#' @return the value of the gradient for the variance function at \code{x} for the difference process \eqn{Z^\Delta = Z_x - \widetilde{Z}_x}.
grad_var_Delta_T<-function(x,kmModel,simupoints,T.mat,F.mat){
  x<-matrix(x,ncol=kmModel@d)
  colnames(x)<-colnames(kmModel@X)
  obj<-krig_weight_GPsimu(object=kmModel,simu_points=simupoints,krig_points=x,T.mat = T.mat,F.mat = F.mat)
  obj_deriv<-grad_kweights(object = kmModel,simu_points = simupoints,krig_points = matrix(x,ncol=kmModel@d),T.mat = T.mat,F.mat = F.mat)
#  krig_mean_init <- matrix(obj_deriv$krig.mean.init,ncol=kmModel@d)
  d_weights <- t(obj_deriv$Lambda.end)
  weights<-t(obj$Lambda.end)
  ggVar<-gradKm_dnewdata(object = kmModel,newdata = x,type = "UK",se.compute = TRUE)$s2

  kEE<-predict.km(object = kmModel,newdata = simupoints,type = "UK",cov.compute = TRUE,light.return = TRUE)$cov

  return(list(ggV=t(ggVar), ggRest= - 2*d_weights%*%tcrossprod(kEE,weights)))
}

# prof_mean_var_Delta function
#' @author Dario Azzimonti
#' @name prof_mean_var_Delta
#' @title Profile extrema for the mean and variance functions of difference process
#' @description The function prof_mean_var_Delta computes the profile extrema functions for the mean and variance functions of the difference process \eqn{Z_x - \widetilde{Z}_x} at \code{x}.
#' @param kmModel the \link[DiceKriging]{km} model of the Gaussian process \eqn{Z}.
#' @param simupoints the matrix \eqn{l x d} containing the simulation points \eqn{E}.
#' @param options_approx an optional list of options for approxMaxMin, see \link{approxMaxMin} for details.
#' @param options_full_sims an optional list of options for getAllMaxMin, see \link{getAllMaxMin} for details. If NULL the full computations are not excuted. NOTE: this computations might be very expensive!
#' @param T.mat the upper triangular factor of the Choleski decomposition of the covariance matrix of \code{rbind(kmModel@X,simupoints)}, if \code{NULL} it is computed.
#' @param F.mat the evaluation of the trend function at \code{rbind(kmModel@X,simupoints)}, see \link[stats]{model.matrix}, if \code{NULL} it is computed.
#' @return the profile extrema functions at \code{options_approx$design} for the mean and variance function of the difference process \eqn{Z^\Delta = Z_x - \widetilde{Z}_x}.
#' @export
prof_mean_var_Delta<-function(kmModel,simupoints,options_full_sims=NULL,options_approx=NULL,F.mat=NULL,T.mat=NULL){
  if(is.null(F.mat))
    F.mat <- model.matrix(object=kmModel@trend.formula, data = data.frame(rbind(kmModel@X,simupoints)))
  if(is.null(T.mat)){
    K <- covMatrix(object=kmModel@covariance,X=rbind(kmModel@X,simupoints))$C
    T.mat <- chol(K)
  }

  ## Set up functions for optimization
  # mean Delta
#  g_mean_spec<-function(x){
#    return(mean_Delta_T(x=x,kmModel=kmModel,simupoints=simupoints,T.mat = T.mat,F.mat = F.mat))
#  }

  # grad mean Delta
#  g_mean_der_spec<-function(x){
#    return(grad_mean_Delta_T(x=x,kmModel=kmModel,simupoints=simupoints,T.mat = T.mat,F.mat = F.mat))
#  }

  # var Delta
  g_var_spec<-function(x){
    temp<-var_Delta_T(x=x,kmModel=kmModel,simupoints=simupoints,T.mat = T.mat,F.mat = F.mat)
    return(temp$vv+temp$vv_E)
  }

  g_var_der_spec<-function(x){
    temp<-grad_var_Delta_T(x=x,kmModel=kmModel,simupoints=simupoints,T.mat = T.mat,F.mat = F.mat)
    return(matrix(temp$ggV+temp$ggRest,ncol=kmModel@d))
  }

  d<-kmModel@d


  if(!is.null(options_full_sims)){
  #  timeIn<-get_nanotime()
  #  mean_delta<-getAllMaxMin(f = g_mean_spec,fprime = g_mean_der_spec,d = d,options = options_full_sims)
  #  tFull_mean<-(get_nanotime()-timeIn)*1e-9

    timeIn<-get_nanotime()
    var_delta<-getAllMaxMin(f = g_var_spec,fprime = g_var_der_spec,d = d,options = options_full_sims)
    time_var<-(get_nanotime()-timeIn)*1e-9

#    times<-c(tFull_mean,tFull_var)

  }else{
 #   timeIn<-get_nanotime()
 #   mean_delta<-approxMaxMin(f = g_mean_spec,fprime = g_mean_der_spec,d = d,opts = options_approx)
 #   tApprox_mean<-(get_nanotime()-timeIn)*1e-9

    timeIn<-get_nanotime()
    var_delta<-approxMaxMin(f = g_var_spec,fprime = g_var_der_spec,d = d,opts = options_approx)
    time_var<-(get_nanotime()-timeIn)*1e-9
  }

  # mean_delta is zero
  mean_delta<-list(res=list(min=matrix(0,ncol=d,nrow=nrow(var_delta$res$min)),
                            max=matrix(0,ncol=d,nrow=nrow(var_delta$res$min))))


  return(list(mean=mean_delta,var=var_delta,times=time_var))
}


# bound_profiles function
#' @author Dario Azzimonti
#' @name bound_profiles
#' @title Bound for profile extrema quantiles
#' @description The function bound_profiles computes the upper and lower bounds for the profile extrema quantiles of a Gaussian process model.
#' @param objectUQ an object returned by \link{coordProf_UQ} or the object saved in \code{obj$res_UQ}, if \code{obj} is the object returned by \link{coordinateProfiles}
#' @param mean_var_delta  the profile extrema functions at \code{options_approx$design} for the mean and variance function of the difference process \eqn{Z^\Delta = Z_x - \widetilde{Z}_x}. Object returned by \link{prof_mean_var_Delta}.
#' @param beta the level of confidence for the approximate simulations
#' @param alpha the level of confidence for the bound
#' @param options_approx an optional list of options for approxMaxMin, see \link{approxMaxMin} for details.
#' @param options_full_sims an optional list of options for getAllMaxMin, see \link{getAllMaxMin} for details. If NULL the full computations are not excuted. NOTE: this computations might be very expensive!
#' @return a list containing \itemize{
#' \item{\code{bound: }}{a list containing the upper/lower bound for profile sup and inf}
#' \item{\code{approx: }}{a list containing the upper/lower approximate quantiles for profile sup and inf}
#' }
#' @export
bound_profiles<-function(objectUQ,mean_var_delta=NULL,beta=0.1,alpha=0.05,options_approx=NULL,options_full_sims=NULL){

  if(is.null(mean_var_delta)){
    mean_var_delta<-prof_mean_var_Delta(kmModel=objectUQ$kmModel,simupoints = objectUQ$sPts$par,options_full_sims=options_full_sims,options_approx=options_approx,F.mat=NULL,T.mat=NULL)
  }

  d<-objectUQ$kmModel@d

  # get (empirical) quantiles from approximations
  approx_quant<-list(lower=list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d),
                              max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d))),
                     upper=list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d),
                                         max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d))) )

  bound_quant<-list(lower=list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d),
                                         max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d))),
                     upper=list(res=list(min=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d),
                                         max=matrix(NA,nrow = options_approx$fullDesignSize,ncol = d))) )

  correction<-log(2/(beta-alpha))

  for(coord in seq(d)){
    approx_quant$lower$res$max[,coord]<-apply(objectUQ$profSups[coord,,],1,function(x){return(quantile(x,beta))})
    bound_quant$lower$res$max[,coord]<-approx_quant$lower$res$max[,coord]+mean_var_delta$mean$res$max[,coord]-sqrt(2*mean_var_delta$var$res$max[,coord]*correction)

    approx_quant$lower$res$min[,coord]<-apply(objectUQ$profInfs[coord,,],1,function(x){return(quantile(x,beta))})
    bound_quant$lower$res$min[,coord]<-approx_quant$lower$res$min[,coord]+mean_var_delta$mean$res$min[,coord]-sqrt(2*mean_var_delta$var$res$min[,coord]*correction)

    approx_quant$upper$res$max[,coord]<-apply(objectUQ$profSups[coord,,],1,function(x){return(quantile(x,1-beta))})
    bound_quant$upper$res$max[,coord]<-approx_quant$upper$res$max[,coord]+mean_var_delta$mean$res$max[,coord]+sqrt(2*mean_var_delta$var$res$max[,coord]*correction)

    approx_quant$upper$res$min[,coord]<-apply(objectUQ$profInfs[coord,,],1,function(x){return(quantile(x,1-beta))})
    bound_quant$upper$res$min[,coord]<-approx_quant$upper$res$min[,coord]+mean_var_delta$mean$res$min[,coord]+sqrt(2*mean_var_delta$var$res$min[,coord]*correction)
  }

  return(list(bound=bound_quant,approx=approx_quant,mean_var_D=mean_var_delta))


}

