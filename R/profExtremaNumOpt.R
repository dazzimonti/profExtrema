#' @author Dario Azzimonti
#' @name getProfileExtrema
#' @title Profile extrema with BFGS optimization
#' @description Evaluate profile extrema for a set of matrices allPsi with full optimization.
#' @param f the function to be evaluated
#' @param fprime derivative of the function
#' @param d dimension of the input domain
#' @param allPsi a list containing the matrices Psi (dim \eqn{pxd}) for which to compute the profile extrema
#' @param opts a list containing the options for this function and the subfunctions \link{getProfileSup_optim}, \link{getProfileInf_optim}. The options only for getProfileExtrema are
#' \itemize{
#' \item{\code{limits:}}{an optional list containing \code{lower} and \code{upper}, two vectors with the limits of the input space. If NULL then \code{limits=list(upper=rep(1,d),lower=rep(0,d))}}
#' \item{\code{discretization:}}{an optional integer representing the discretization size for the profile computation for each dimension of eta. Pay attention that this leads to a grid of size \code{discretization^p}.}
#' \item{\code{heavyReturn:}}{If TRUE returns also all minimizers, default is FALSE.}
#' \item{\code{plts:}}{If TRUE and p==1 for all Psi in allPsi, plots the profile functions at each Psi, default is FALSE.}
#' \item{\code{verb:}}{If TRUE, outputs intermediate results, default is FALSE.}
#' }
#' @return a list of two data frames (min, max) of the evaluations of \eqn{P^sup_Psi f(eta) = sup_{Psi x = \eta} f(x) } and \eqn{P^inf_Psi f(eta) = inf_{Psi x = \eta} f(x) }
#' discretized over 50 equally spaced points for each dimension for each Psi in \code{allPsi}. This number can be changed by defining it in options$discretization.
#' @examples
#' # Compute the oblique profile extrema with full optimization on 2d example
#'
#' # Define the function
#' testF <- function(x,params,v1=c(1,0),v2=c(0,1)){
#' return(sin(crossprod(v1,x)*params[1]+params[2])+cos(crossprod(v2,x)*params[3]+params[4])-1.5)
#' }
#'
#' testFprime <- function(x,params,v1=c(1,0),v2=c(0,1)){
#'   return(matrix(c(params[1]*v1[1]*cos(crossprod(v1,x)*params[1]+params[2])-
#'                   params[3]*v2[1]*sin(crossprod(v2,x)*params[3]+params[4]),
#'                  params[1]*v1[2]*cos(crossprod(v1,x)*params[1]+params[2])-
#'                   params[3]*v2[2]*sin(crossprod(v2,x)*params[3]+params[4])),ncol=1))
#' }
#'
#'
#' # Define the main directions of the function
#' theta=pi/6
#' pparams<-c(1,0,10,0)
#' vv1<-c(cos(theta),sin(theta))
#' vv2<-c(cos(theta+pi/2),sin(theta+pi/2))
#'
#' # Define optimizer friendly function
#' f <-function(x){
#' return(testF(x,pparams,vv1,vv2))
#' }
#' fprime <- function(x){
#'  return(testFprime(x,pparams,vv1,vv2))
#' }
#'
#' # Define list of directions where to evaluate the profile extrema
#' all_Psi <- list(Psi1=vv1,Psi2=vv2)
#'
#'
#' \donttest{
#' # Evaluate profile extrema along directions of all_Psi
#' allOblique<-getProfileExtrema(f=f,fprime = fprime,d = 2,allPsi = all_Psi,
#'                               opts = list(plts=FALSE,discretization=100,multistart=8))
#'
#'
#' # Consider threshold=0
#' threshold <- 0
#'
#' # Plot oblique profile extrema functions
#' plotMaxMin(allOblique,allOblique$Design,threshold = threshold)
#'
#' ## Since the example is two dimensional we can visualize the regions excluded by the profile extrema
#' # evaluate the function at a grid for plots
#' inDes<-seq(0,1,,100)
#' inputs<-expand.grid(inDes,inDes)
#' outs<-apply(X = inputs,MARGIN = 1,function(x){return(testF(x,pparams,v1=vv1,v2=vv2))})
#'
#' # obtain the points where the profiles take the threshold value
#' cccObl<-getChangePoints(threshold = threshold,allRes = allOblique,Design = allOblique$Design)
#'
#' # visualize the functions and the regions excluded
#'
#' image(inDes,inDes,matrix(outs,ncol=100),col=grey.colors(20),main="Example and oblique profiles")
#' contour(inDes,inDes,matrix(outs,ncol=100),add=T,nlevels = 20)
#' contour(inDes,inDes,matrix(outs,ncol=100),add=T,levels = c(threshold),col=4,lwd=1.5)
#' plotOblique(cccObl$alwaysEx$`0`[[1]],all_Psi[[1]],col=3)
#' plotOblique(cccObl$alwaysEx$`0`[[2]],all_Psi[[2]],col=3)
#' plotOblique(cccObl$neverEx$`0`[[1]],all_Psi[[1]],col=2)
#' plotOblique(cccObl$neverEx$`0`[[2]],all_Psi[[2]],col=2)
#'
#' }
#' @export
getProfileExtrema<-function(f,fprime=NULL,d,allPsi,opts=NULL){

  num_Psi<-length(allPsi)

  # initialize opts
  if(is.null(opts$limits)){
    ll_b = rep(0,d)
    uu_b = rep(1,d)
  }else{
    ll_b = opts$limits$lower
    uu_b = opts$limits$upper
  }

  if(is.null(opts$discretization)){
    dd_eta <- 50
  }else{
    dd_eta <- opts$discretization
  }

  if(is.null(opts$multistart)){
    opts$multistart=5
  }

  # Useful for choosing limits for eta
  cubeVertex<-matrix(c(ll_b[1],uu_b[1]),nrow=1)
  ii=1
  while(ii<d){
    ii=ii+1
    cubeVertex <- cbind(cubeVertex,cubeVertex)
    cubeVertex <-rbind(cubeVertex, c(rep(ll_b[ii],ncol(cubeVertex)/2),c(rep(uu_b[ii],ncol(cubeVertex)/2))))
  }

  p=nrow(matrix(allPsi[[1]],ncol=d))

  allMaxPoints<-list()
  allMinPoints<-list()
  results<-list(min=data.frame(matrix(NA,ncol=num_Psi,nrow = dd_eta^p)),max=data.frame(matrix(NA,ncol=num_Psi,nrow = dd_eta^p)))

  # Save the design, useful for plotting functions
  Design<- matrix(NA,ncol=num_Psi,nrow=dd_eta*p)

  # Loop over the different Psi
  for(i in seq(num_Psi)){

    if(!is.null(opts$verb))
      if(opts$verb)
        cat("Psi number ",i," of ",num_Psi,"\n")

    # take current Psi
    cPsi = allPsi[[i]]

    p = nrow(cPsi)
    if(is.null(p)){
      cPsi <-matrix(cPsi,ncol=d)
      p=nrow(cPsi)
    }
    # Choose limits for etas for current Psi
    mmEtas<-apply(crossprod(t(cPsi),cubeVertex),1,min)
    MMetas<-apply(crossprod(t(cPsi),cubeVertex),1,max)

    if(p==1){
      etas1<-matrix(seq(from = mmEtas,to = MMetas,length.out = dd_eta),ncol=1)
      Design[,i]<-etas1
    }else if(p==2){
      etas1<-expand.grid(seq(from = mmEtas[1],to = MMetas[1],length.out = dd_eta),seq(from = mmEtas[2],to = MMetas[2],length.out = dd_eta))
      Design[,i]<-c(seq(from = mmEtas[1],to = MMetas[1],length.out = dd_eta),seq(from = mmEtas[2],to = MMetas[2],length.out = dd_eta))
    }else{
      stop("ERROR: no implementation for p>2!")
    }

    pSup<-apply(etas1,1,function(x){a_res=getProfileSup_optim(eta = x,Psi = matrix(cPsi,ncol=d),f = f,fprime = fprime,d = d,options = opts);gc();return(a_res)})
    pInf<-apply(etas1,1,function(x){a_res=getProfileInf_optim(eta = x,Psi = matrix(cPsi,ncol=d),f = f,fprime = fprime,d = d,options = opts);gc();return(a_res)})

    results$max[[i]] <- sapply(pSup,function(x){x$val})
    allMaxPoints[[i]]<- sapply(pSup,function(x){x$aux$solution})

    results$min[[i]] <- sapply(pInf,function(x){x$val})
    allMinPoints[[i]]<- sapply(pSup,function(x){x$aux$solution})



    if(!is.null(opts$plts)){
      if(opts$plts){
        if(p==1){
          plot(etas1,results$min[,i],ylim=c(min(results$min[,i]),max(results$max[,i])),type='l',main=paste("Psi number",i))
          lines(etas1,results$max[,i])
        }else if(p==2){
          par(mfrow=c(1,2))
          image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(results$max[[i]],ncol=dd_eta),col=gray.colors(20),main=bquote(P[Psi[.(i)]]^sup ~ "f"))
          contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(results$max[[i]],ncol=dd_eta),nlevels = 10,add=T)

          image(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(results$min[[i]],ncol=dd_eta),col=gray.colors(20),main=bquote(P[Psi[.(i)]]^inf ~ "f"))
          contour(x=Design[1:(dd_eta),i],y=Design[(dd_eta+1):(2*dd_eta),i],z=matrix(results$min[[i]],ncol=dd_eta),nlevels = 10,add=T)
          par(mfrow=c(1,1))
        }else{
          stop("We should never be here!")
        }
      }
    }


  }

  toReturn<-list(res=results)
  toReturn$Design<-Design
  if(is.null(opts$heavyReturn))
    opts$heavyReturn=FALSE
  if(opts$heavyReturn){
    toReturn$minima<-allMinPoints
    toReturn$maxima<-allMaxPoints
  }
  return(toReturn)
}


# approxProfileExtrema
#' @author Dario Azzimonti
#' @name approxProfileExtrema
#' @title Approximate profile extrema functions
#' @description Evaluate profile extrema for a set of Psi with approximations at few values
#' @param f the function to be evaluated
#' @param fprime derivative of the function
#' @param d dimension of the input domain
#' @param allPsi a list containing the matrices Psi (dim \eqn{pxd}) for which to compute the profile extrema
#' @param opts a list containing the options for this function and the subfunctions \link{getProfileSup_optim}, \link{getProfileInf_optim} or \link{getProfileExtrema}. The options only for approxProfileExtrema are
#' \itemize{
#' \item{\code{limits:}}{an optional list with the upper and lower limits of input space dimension, if NULL then \code{limits=list(upper=rep(1,d),lower=rep(0,d))}}
#' \item{\code{smoother:}}{Select which smoother to use:a string that selects which smoother to use: \itemize{
#'       \item{\code{"1order"}}: first order interpolation with gradient
#'       \item{\code{"splineSmooth"}}: smoothing spline with default degrees of freedom (DEFAULT OPTION)
#'       \item{\code{"quantSpline"}}: profile inf and profile sup approximated with quantile spline regression at levels 0.1 and 0.9 respectively
#' }}
#' \item{\code{heavyReturn:}}{If TRUE returns also all minimizers, default is FALSE.}
#' \item{\code{initDesign:}}{A list of the same length as allPsi containing the designs of few points where the expensive sup is evaluated. If Null it is automatically initialized}
#' \item{\code{fullDesignSize:}}{The full design where the function is approximated.}
#' \item{\code{multistart:}}{number of multistarts for optim procedure.}
#' \item{\code{numMCsamples:}}{number of MC samples for the sup.}
#' \item{\code{plts:}}{If TRUE, plots the max/min functions at each coordinate, default is FALSE.}
#' \item{\code{verb:}}{If TRUE, outputs intermediate results, default is FALSE.}
#' }
#' @return a list of two data frames (min, max) of the evaluations of \eqn{f_sup(x_i) = sup_{x_j \neq i} f(x_1,\dots,x_d) } and \eqn{f_inf(x_i) = inf_{x_j \neq i} f(x_1,\dots,x_d) }
#' for each i at the design Design. By default Design is a 100 equally spaced points for each dimension. It can be changed by defining it in options$Design
#' @examples
#' # Compute the oblique profile extrema with approximate optimization on 2d example
#'
#' # Define the function
#' testF <- function(x,params,v1=c(1,0),v2=c(0,1)){
#' return(sin(crossprod(v1,x)*params[1]+params[2])+cos(crossprod(v2,x)*params[3]+params[4])-1.5)
#' }
#'
#' testFprime <- function(x,params,v1=c(1,0),v2=c(0,1)){
#'   return(matrix(c(params[1]*v1[1]*cos(crossprod(v1,x)*params[1]+params[2])-
#'                   params[3]*v2[1]*sin(crossprod(v2,x)*params[3]+params[4]),
#'                  params[1]*v1[2]*cos(crossprod(v1,x)*params[1]+params[2])-
#'                   params[3]*v2[2]*sin(crossprod(v2,x)*params[3]+params[4])),ncol=1))
#' }
#'
#'
#' # Define the main directions of the function
#' theta=pi/6
#' pparams<-c(1,0,10,0)
#' vv1<-c(cos(theta),sin(theta))
#' vv2<-c(cos(theta+pi/2),sin(theta+pi/2))
#'
#' # Define optimizer friendly function
#' f <-function(x){
#' return(testF(x,pparams,vv1,vv2))
#' }
#' fprime <- function(x){
#'  return(testFprime(x,pparams,vv1,vv2))
#' }
#'
#' # Define list of directions where to evaluate the profile extrema
#' all_Psi <- list(Psi1=vv1,Psi2=vv2)
#'
#' # Evaluate profile extrema along directions of all_Psi
#' allOblique<-approxProfileExtrema(f=f,fprime = fprime,d = 2,allPsi = all_Psi,
#'                                  opts = list(plts=FALSE,heavyReturn=TRUE))
#'
#' \donttest{
#' # Consider threshold=0
#' threshold <- 0
#'
#' # Plot oblique profile extrema functions
#' plotMaxMin(allOblique,allOblique$Design,threshold = threshold)
#'
#' ## Since the example is two dimensional we can visualize the regions excluded by the profile extrema
#' # evaluate the function at a grid for plots
#' inDes<-seq(0,1,,100)
#' inputs<-expand.grid(inDes,inDes)
#' outs<-apply(X = inputs,MARGIN = 1,function(x){return(testF(x,pparams,v1=vv1,v2=vv2))})
#'
#' # obtain the points where the profiles take the threshold value
#' cccObl<-getChangePoints(threshold = threshold,allRes = allOblique,Design = allOblique$Design)
#'
#' # visualize the functions and the regions excluded
#'
#' image(inDes,inDes,matrix(outs,ncol=100),col=grey.colors(20),main="Example and oblique profiles")
#' contour(inDes,inDes,matrix(outs,ncol=100),add=T,nlevels = 20)
#' contour(inDes,inDes,matrix(outs,ncol=100),add=T,levels = c(threshold),col=4,lwd=1.5)
#' plotOblique(cccObl$alwaysEx$`0`[[1]],all_Psi[[1]],col=3)
#' plotOblique(cccObl$alwaysEx$`0`[[2]],all_Psi[[2]],col=3)
#' plotOblique(cccObl$neverEx$`0`[[1]],all_Psi[[1]],col=2)
#' plotOblique(cccObl$neverEx$`0`[[2]],all_Psi[[2]],col=2)
#'
#' }
#' @export
approxProfileExtrema=function(f,fprime=NULL,d,allPsi,opts=NULL){

  # Set up options
  if(is.null(opts$limits)){
    limits<-list(lower=rep(0,d),upper=rep(1,d))
  }

  num_Psi<-length(allPsi)

  # initialize opts
  if(is.null(opts$limits)){
    ll_b = rep(0,d)
    uu_b = rep(1,d)
  }else{
    ll_b = limits$lower
    uu_b = limits$upper
  }

  if(is.null(opts$fullDesignSize)){
    opts$fullDesignSize<-80*d
  }

  if(is.null(opts$multistart)){
    opts$multistart=5
  }
  if(is.null(opts$initDesign))
    opts$initDesign<-list()

  # Useful for choosing limits for eta
  cubeVertex<-matrix(c(ll_b[1],uu_b[1]),nrow=1)
  ii=1
  while(ii<d){
    ii=ii+1
    cubeVertex <- cbind(cubeVertex,cubeVertex)
    cubeVertex <-rbind(cubeVertex, c(rep(ll_b[ii],ncol(cubeVertex)/2),c(rep(uu_b[ii],ncol(cubeVertex)/2))))
  }

  p=nrow(matrix(allPsi[[1]],ncol=d))

  exact_allMaxPoints<-list()
  exact_allMinPoints<-list()
  exact_results<-list(min=data.frame(matrix(NA,ncol=num_Psi,nrow = ceiling(sqrt(d*p)*10))),max=data.frame(matrix(NA,ncol=num_Psi,nrow = ceiling(sqrt(d*p)*10))))

  # initialize the max/min values
  sK_max<-data.frame(matrix(NA,nrow=opts$fullDesignSize^p,ncol=num_Psi))
  sK_min<-data.frame(matrix(NA,nrow=opts$fullDesignSize^p,ncol=num_Psi))

  # Save Design to return
  Design<- matrix(NA,ncol=num_Psi,nrow=opts$fullDesignSize*p)

  # Loop over the different Psi
  for(i in seq(num_Psi)){

    if(!is.null(opts$verb))
      if(opts$verb)
        cat("Psi number ",i," of ",num_Psi,"\n")

    # take current Psi
    cPsi = allPsi[[i]]

    p = nrow(cPsi)
    if(is.null(p)){
      cPsi <-matrix(cPsi,ncol=d)
      p=nrow(cPsi)
    }
    # Choose limits for etas for current Psi
    mmEtas<-apply(crossprod(t(cPsi),cubeVertex),1,min)
    MMetas<-apply(crossprod(t(cPsi),cubeVertex),1,max)


    # Get initial design
    if(length(opts$initDesign)<num_Psi){
      if(p==1){
        opts$initDesign[[i]]<-matrix(seq(from=mmEtas[1],to=MMetas[1],,ceiling(sqrt(d)*10)),ncol=1)
      }else{
        ### NEEDS TEST!
        opts$initDesign[[i]]<-matrix(mmEtas,ncol=p,byrow = T,nrow=ceiling(sqrt(d*p)*10))+maximinLHS(ceiling(sqrt(d*p)*10),p)%*%diag(MMetas-mmEtas,ncol=p)
      }

    }


    pSup<-apply(opts$initDesign[[i]],1,function(x){return(getProfileSup_optim(eta = x,Psi = matrix(cPsi,ncol=d),f = f,fprime = fprime,d = d,options = opts))})
    gc()
    pInf<-apply(opts$initDesign[[i]],1,function(x){return(getProfileInf_optim(eta = x,Psi = matrix(cPsi,ncol=d),f = f,fprime = fprime,d = d,options = opts))})
    gc()

    exact_results$max[[i]] <- sapply(pSup,function(x){x$val})
    exact_allMaxPoints[[i]]<- sapply(pSup,function(x){x$aux$solution})
    exact_results$min[[i]] <- sapply(pInf,function(x){x$val})
    exact_allMinPoints[[i]]<- sapply(pInf,function(x){x$aux$solution})

    # the gradient of the sup/inf function is saved
    if(!is.null(fprime) & p==1){
      nn<-length(exact_results$max[[i]])
      exact_grad_min<-matrix(NA,ncol=1,nrow=nn)
      exact_grad_max<-matrix(NA,ncol=1,nrow=nn)
      for(ll in seq(nn)){

        yy<-exact_allMinPoints[[i]][,ll]
        exact_grad_min[ll,]<-fprime(yy)[i,] ### THIS IS WRONG!
        yy<-exact_allMaxPoints[[i]][,ll]
        exact_grad_max[ll,]<-fprime(yy)[i,] ### THIS IS WRONG!
      }
    }
    if(p==1){
      newPoints<-matrix(seq(mmEtas,MMetas,,opts$fullDesignSize),ncol=1)
      Design[,i]<-newPoints
    }else{
      newPoints<-expand.grid(seq(mmEtas[1],MMetas[1],,opts$fullDesignSize),seq(mmEtas[2],MMetas[2],,opts$fullDesignSize))
      Design[,i]<-c(seq(mmEtas[1],MMetas[1],,opts$fullDesignSize),seq(mmEtas[2],MMetas[2],,opts$fullDesignSize))
    }


    if(p==1){

      if(!is.null(opts$smoother) && opts$smoother=="1order"){
        #    sK_max[,i]<- kGradSmooth(newPoints=newPoints,profPoints = exact_allMaxPoints[[i]],
        #                                 profEvals=exact_results$max[[i]], profGradient=exact_grad_max)

        #      sK_min[,i]<-kGradSmooth(newPoints=newPoints,profPoints = exact_allMaxPoints[[i]],
        #                                 profEvals=exact_results$min[[i]], profGradient=exact_grad_min)
        warning(paste("The method ",opts$smoother," is not implemented yet!"))
      }else if(!is.null(opts$smoother) && opts$smoother=="quantSpline"){
        dd_fit<-data.frame(x=opts$initDesign[[i]],y=exact_results$max[[i]])
        fit.max <- rq(y ~ bs(x, df=max(nn%/%1.5,3)),data=dd_fit,tau=0.9)
        dd_fit<-data.frame(x=opts$initDesign[[i]],y=exact_results$min[[i]])
        fit.min <- rq(y ~ bs(x, df=max(nn%/%1.5,3)), tau=0.1,data=dd_fit)
        newPoints<-data.frame(x=newPoints)

        sK_max[,i]<- predict.rq(object = fit.max,newdata = newPoints)
        sK_min[,i]<- predict.rq(object = fit.min,newdata = newPoints)
      }else {
        sk_maxSmooth<-smooth.spline(x=opts$initDesign[[i]],y=exact_results$max[[i]])
        sk_minSmooth<-smooth.spline(x=opts$initDesign[[i]],y=exact_results$min[[i]])

        sK_max[,i]<-predict(sk_maxSmooth,newPoints)$y
        sK_min[,i]<-predict(sk_minSmooth,newPoints)$y
      }

    }else{
      approxKm<-km(design = opts$initDesign[[i]],response = exact_results$max[[i]],covtype = "matern3_2",control = list(trace=F))
      sK_max[,i]<-predict.km(object = approxKm,newdata = newPoints,type = "UK",checkNames = F)$mean

      approxKm<-km(design = opts$initDesign[[i]],response = exact_results$min[[i]],covtype = "matern3_2",control = list(trace=F))
      sK_min[,i]<-predict.km(object = approxKm,newdata = newPoints,type = "UK",checkNames = F)$mean
    }


    ## Debug plots
    if(!is.null(opts$debug) && opts$debug){

      warning("The debug plots do not work!")

      # par(mfrow=c(1,2))
      # image(x=Design[1:opts$fullDesignSize,i],y=Design[(opts$fullDesignSize+1):(opts$fullDesignSize*2),i],
      #       z=matrix(sK_max[[i]],ncol=opts$fullDesignSize),col=gray.colors(20),main=sprintf("Psup Psi %i",i))
      # contour(x=Design[1:opts$fullDesignSize,i],y=Design[(opts$fullDesignSize+1):(opts$fullDesignSize*2),i],
      #       z=matrix(sK_max[[i]],ncol=opts$fullDesignSize),nlevels = 10,add=T)
      # contour(x=Design[1:opts$fullDesignSize,i],y=Design[(opts$fullDesignSize+1):(opts$fullDesignSize*2),i],
      #         z=matrix(sK_max[[i]],ncol=opts$fullDesignSize),levels = threshold,add=T,col=2)
      #
      # image(x=Design[1:opts$fullDesignSize,i],y=Design[(opts$fullDesignSize+1):(opts$fullDesignSize*2),i],
      #       z=matrix(sK_min[[i]],ncol=opts$fullDesignSize),col=gray.colors(20),main=sprintf("Pinf Psi %i",i))
      # contour(x=Design[1:opts$fullDesignSize,i],y=Design[(opts$fullDesignSize+1):(opts$fullDesignSize*2),i],
      #         z=matrix(sK_min[[i]],ncol=opts$fullDesignSize),nlevels = 10,add=T)
      # contour(x=Design[1:opts$fullDesignSize,i],y=Design[(opts$fullDesignSize+1):(opts$fullDesignSize*2),i],
      #         z=matrix(sK_min[[i]],ncol=opts$fullDesignSize),levels = threshold,add=T,col=2)
      # par(mfrow=c(1,1))
      #    par(mfrow=c(1,2))
      # Plot  max/min function
      # ylimTemp<-range(c(sK_min[,coord],aa$res$min[,coord],sK_max[,coord],aa$res$max[,coord]))
      # plot(aa$minima[((coord-1)*nn+1):(coord*nn),coord],aa$res$min[,coord],main=paste("Coordinate",coord,", Min"),
      #      xlab="x",ylab="f",ylim=ylimTemp)
      # lines(newPoints,sK_min[,coord],type='l',col=3)
      # points(aa$maxima[((coord-1)*nn+1):(coord*nn),coord],aa$res$max[,coord])
      # lines(newPoints,sK_max[,coord],type='l',col=3)
      #
      # if(!is.null(opts$smooth) && opts$smooth){
      #   points(sk_minSmooth$x,sk_minSmooth$y,col=2)
      #   points(sk_maxSmooth$x,sk_maxSmooth$y,col=2)
      # }
    }
  }

  res=list(min=sK_min,max=sK_max)
  toReturn<-list(res=res)
  toReturn$Design <- Design
  if(!is.null(opts$heavyReturn) && opts$heavyReturn){
    aa<-list(e_res=exact_results,e_maxPts= exact_allMaxPoints,e_minPts=exact_allMinPoints)
    aa$design<-opts$initDesign
    toReturn$profPoints<-aa
  }
  return(toReturn)



}

