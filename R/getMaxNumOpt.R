# getMax function
#' @author Dario Azzimonti
#' @name getMax
#' @title Coordinate profile sup function
#' @description Compute coordinate profile sup functions
#' @param x one dimensional point where the function is to be evaluated
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format)
#' @param coord integer selecting the dimension that is fixed, the other ones are optimized over
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to optim:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d-1)}
#' \item{\code{method:}}{ is the string denoting the chosen method for the optimization (see optim for details)}
#' \item{\code{lower:} }{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{max_{x_1,\dots, x_{coord-1},x_{coord+1}, \dots, x_d} f(x_1,\dots,x_d)}
getMax = function(x,f,fprime,coord,d,options=NULL){
  fxFixed=function(y){
    res<-NA
    if(coord==1)
      res<-f(c(x,y))
    if(coord==d)
      res<-f(c(y,x))
    if(coord>1 && coord<d)
      res<-f(c(y[1:(coord-1)],x,y[coord:(d-1)]))
    return(res)
  }
  if(!is.null(fprime)){
    fprimexFixed=function(y){
      res<-NA
      if(coord==1)
        res<-fprime(c(x,y))
      if(coord==d)
        res<-fprime(c(y,x))
      if(coord>1 && coord<d)
        res<-fprime(c(y[1:(coord-1)],x,y[coord:(d-1)]))
      return(res[-coord])
    }
  }else{
    fprimexFixed<-NULL
  }

  startingPoint<-runif(d-1)
  if(!is.null(options$par)){
    startingPoint<-options$par
  }

#  if(d<=2){
#    methd="Brent"
#    lower=0
#    upper=1
#  }else{
    methd="L-BFGS-B"
    lower=rep(0,d-1)
    upper=rep(1,d-1)
#  }
  trace=0

  if(!is.null(options$method))
    methd=options$method
  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace


  optRes<-list(value=-Inf,aux=NA)
  if(!is.null(options$multistart)){
    optRess<-NULL
    for(i in seq(options$multistart)){
      startingPoint<-runif(d-1,min = lower,max=upper)
      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
      if(optRess$value>=optRes$value){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    optRes<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
  }



  return(list(val=optRes$value,aux=optRes))

}

# getMin function
#' @author Dario Azzimonti
#' @name getMin
#' @title Coordinate profile inf function
#' @description Compute coordinate profile inf functions
#' @param x one dimensional point where the function is to be evaluated
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format)
#' @param coord integer selecting the dimension that is fixed, the other ones are optimized over
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to optim:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d-1)}
#' \item{\code{method:}}{is the string denoting the chosen method for the optimization (see optim for details)}
#' \item{\code{lower:}}{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{min_{x_1,\dots, x_{coord-1},x_{coord+1}, \dots, x_d} f(x_1,\dots,x_d)}
getMin = function(x,f,fprime,coord,d,options=NULL){
  fxFixed=function(y){
    res<-NA
    if(coord==1)
      res<-f(c(x,y))
    if(coord==d)
      res<-f(c(y,x))
    if(coord>1 && coord<d)
      res<-f(c(y[1:(coord-1)],x,y[coord:(d-1)]))
    return(res)
  }
  if(!is.null(fprime)){
    fprimexFixed=function(y){
      res<-NA
      if(coord==1)
        res<-fprime(c(x,y))
      if(coord==d)
        res<-fprime(c(y,x))
      if(coord>1 && coord<d)
        res<-fprime(c(y[1:(coord-1)],x,y[coord:(d-1)]))
      return(res[-coord])
    }
  }else{
    fprimexFixed<-NULL
  }

  startingPoint<-runif(d-1)
  if(!is.null(options$par)){
    startingPoint<-options$par
  }

  # if(d<=2){
  #   methd="Brent"
  #   lower=0
  #   upper=1
  # }else{
  methd="L-BFGS-B"
  lower<-rep(0,d-1)
  upper<-rep(1,d-1)
  # }
  trace=0

  if(!is.null(options$method))
    methd=options$method
  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace

  optRes<-list(value=Inf,aux=NA)
  if(!is.null(options$multistart)){
    optRess<-NULL
    for(i in seq(options$multistart)){
      startingPoint<-runif(d-1,min = lower,max=upper)
      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(trace=trace),method = methd,lower = lower,upper = upper)
      if(optRess$value<optRes$value){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    optRes<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(trace=trace),method = methd,lower = lower,upper = upper)
  }

  return(list(val=optRes$value,aux=optRes))

}


# getMaxMinMC function
#' @author Dario Azzimonti
#' @name getMaxMinMC
#' @title Coordinate profile extrema with MC
#' @description Compute coordinate profile extrema with MC
#' @param x one dimensional point where the function is to be evaluated
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format)
#' @param coord integer selecting the dimension that is fixed, the other ones are optimized over
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to the MC optimizer:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d-1)}
#' \item{\code{numMCsamples:}}{ number of MC samples}
#' \item{\code{rand}}{string that chooses the type of randomness in MC: "unif" (uniform in [lower,upper]), "norm" (independent normal with mean 0 and variance 1)}
#' \item{\code{lower:} }{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{max_{x_1,\dots, x_{coord-1},x_{coord+1}, \dots, x_d} f(x_1,\dots,x_d)}
getMaxMinMC = function(x,f,fprime,coord,d,options=NULL){
  fxFixed=function(y){
    res<-NA
    if(coord==1)
      res<-f(c(x,y))
    if(coord==d)
      res<-f(c(y,x))
    if(coord>1 && coord<d)
      res<-f(c(y[1:(coord-1)],x,y[coord:(d-1)]))
    return(res)
  }
  if(!is.null(fprime)){
    fprimexFixed=function(y){
      res<-NA
      if(coord==1)
        res<-fprime(c(x,y))
      if(coord==d)
        res<-fprime(c(y,x))
      if(coord>1 && coord<d)
        res<-fprime(c(y[1:(coord-1)],x,y[coord:(d-1)]))
      return(res[-coord])
    }
  }else{
    fprimexFixed<-NULL
  }

  startingPoint<-runif(d-1)
  if(!is.null(options$par)){
    startingPoint<-options$par
  }


  lower=rep(0,d-1)
  upper=rep(1,d-1)
  numMCsamples=1000
  trace=0

  if(!is.null(options$numMCsamples))
    numMCsamples=options$numMCsamples
  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace


  optRes<-list(value=-Inf,aux=NA)

  # sample numSamplesMC/2 uniformly(/sobol?) in domain
  # evaluate f and fprime there
  # choose
#  for(i in seq(options$)){
  if(is.null(options$rand) | options$rand=="unif"){
    firstPassage<-matrix(runif((d-1)*numMCsamples,min=lower,max=upper),ncol=(d-1))
  }else if(!is.null(options$rand) & options$rand=="norm"){
    firstPassage<-matrix(rnorm((d-1)*numMCsamples),ncol=(d-1))
  }
  firstRes<-apply(firstPassage,1,fxFixed)

  idxMax<-which.max(firstRes)
  idxMin<-which.min(firstRes)


  return(list(max=list(val=firstRes[idxMax],pts=firstPassage[idxMax,]),
              min=list(val=firstRes[idxMin],pts=firstPassage[idxMin,]),
              aux=list(points=firstPassage,vals=firstRes)))

}


#' @author Dario Azzimonti
#' @name getAllMaxMin
#' @title Coordinate profile extrema with BFGS optimization
#' @description Evaluate coordinate profile extrema with full optimization.
#' @param f the function to be evaluated
#' @param fprime derivative of the function
#' @param d dimension of the input domain
#' @param options a list containing the options for this function and the subfunctions getMax, getMin
#' see documentation of getMax, getMin for details. The options only for getAllMaxMin are
#' \itemize{
#' \item{\code{Design:}}{an optional design matrix with the discretization of each dimension, if NULL then for each dimension Design[,coord] = seq(0,1,length.out=100)}
#' \item{\code{heavyReturn:}}{If TRUE returns also all minimizers, default is FALSE.}
#' \item{\code{plts:}}{If TRUE, plots the max/min functions at each coordinate, default is FALSE.}
#' \item{\code{verb:}}{If TRUE, outputs intermediate results, default is FALSE.}
#' \item{\code{MonteCarlo:}}{If TRUE, use the MC optimizer otherwise use standard optim.}
#' }
#' @return a list of two data frames (min, max) of the evaluations of \eqn{f_sup(x_i) = sup_{x_j \neq i} f(x_1,\dots,x_d) } and \eqn{f_inf(x_i) = inf_{x_j \neq i} f(x_1,\dots,x_d) }
#' for each i at the design Design. By default Design is a 100 equally spaced points for each dimension. It can be changed by defining it in options$Design
#' @examples
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # Compute the coordinate profile extrema with full optimization on 2d example
#'
#' # Define the function
#' g=function(x){
#'   return(-branin(x))
#' }
#' # Define the gradient
#' gprime = function(x){
#'   x1 = x[1]*15-5
#'   x2 = x[2]*15
#'   f1prime = (15*25)/(4*pi^4)*x1^3 - (15*75)/(2*pi^3)*x1^2 +
#'   (80*15)/(pi^2)*x1 - (5*15)/(pi^2)*x2*x1 +
#'   10*15/pi*x2 - 60*15/pi-10*15* (1 - 1/(8*pi))*sin(x1)
#'   f2prime = 2*15*(x2-5/(4*pi^2)*x1^2 +5/pi*x1-6)
#'   return(c(-f1prime,-f2prime))
#' }
#' # set up dimension
#' coordProf<-getAllMaxMin(f = g,fprime = gprime,d=2,options = list(multistart=4,heavyReturn=TRUE))
#'
#' \donttest{
#' # Consider threshold=-10
#' threshold<- -10
#' # obtain the points where the profiles take the threshold value
#' pp_change<-getChangePoints(threshold = threshold,allRes = coordProf)
#' # evaluate g at a grid and plot the image
#' x<-seq(0,1,,100)
#' grid<-expand.grid(x,x)
#' g_evals<- apply(X = grid,MARGIN = 1,FUN = g)
#' image(x = x,y = x,z = matrix(g_evals,nrow = 100),col = grey.colors(20))
#' contour(x=x,y=x,z=matrix(g_evals,nrow = 100), add=TRUE, nlevels = 20)
#' contour(x=x,y=x,z=matrix(g_evals,nrow = 100), add=TRUE, levels = threshold,col=2)
#' abline(h = pp_change$neverEx$`-10`[[2]],col="darkgreen",lwd=2)
#' abline(v = pp_change$neverEx$`-10`[[1]],col="darkgreen",lwd=2)
#' # Plot the coordinate profiles and a threshold
#' plotMaxMin(allRes = coordProf,threshold = threshold,changes = TRUE)
#' }
#' @export
getAllMaxMin<-function(f,fprime=NULL,d,options=NULL){

  if(is.null(options$Design)){
    warning("Design not chosen, by default [0,1]^d")
    Design<-replicate(d,seq(0,1,,100))
  }else{
    Design<-options$Design
  }
  nPtsPerDim<-nrow(Design)

  allMaxPoints<-matrix(NA,nrow=(nPtsPerDim*d),ncol=d)
  allMinPoints<-matrix(NA,nrow=(nPtsPerDim*d),ncol=d)
  results<-list(min=data.frame(matrix(NA,ncol=d,nrow = nPtsPerDim)),max=data.frame(matrix(NA,ncol=d,nrow = nPtsPerDim)))
  normGradRes<-array(NA,dim = c(nPtsPerDim,d,2))
  gradOpt<-array(NA,dim = c(nPtsPerDim,d,2))
#  tempMCvals<-list()
#  tempMCpoints<-list()
  for(coord in seq(d)){
    if(!is.null(options$verb))
      if(options$verb)
        cat("Coordinate ",coord," of ",d,"\n")
#    tempMCvals[[coord]]<-matrix(NA,ncol=nPtsPerDim,nrow=options$numMCsamples)
#    tempMCpoints[[coord]]<-array(NA,dim = c(nPtsPerDim,options$numMCsamples,d))
    for(i in seq(nPtsPerDim)){
      if(!is.null(options$verb))
        if(options$verb)
          if(i%%10==0)
            cat("Iteration ",i," of ",nPtsPerDim,"\n")

      if(!is.null(options$MonteCarlo) && options$MonteCarlo){
        temp<-getMaxMinMC(x=Design[i,coord],f = f,fprime = fprime,coord,d,options=options)
#        tempMCvals[[coord]][,i]<-temp$aux$vals
#        pp<-matrix(NA,ncol=d,nrow=options$numMCsamples)
#        pp[,-coord]<-temp$aux$points
#        pp[,coord]<-rep(Design[i,coord],options$numMCsamples)
#        tempMCpoints[[coord]][i,,]<-pp
        results$max[i,coord]<-temp$max$val
        allMaxPoints[(coord-1)*nPtsPerDim+i,-coord]<-temp$max$pts
        allMaxPoints[(coord-1)*nPtsPerDim+i,coord]<-Design[i,coord]

        results$min[i,coord]<-temp$min$val
        allMinPoints[(coord-1)*nPtsPerDim+i,-coord]<-temp$min$pts
        allMinPoints[(coord-1)*nPtsPerDim+i,coord]<-Design[i,coord]

      }else{
        options$par<-allMaxPoints[getClosePoints(x = Design[i,coord],allPts = allMaxPoints,whichDim = coord),-coord]
        options$par<-options$par*(1+rnorm(n=1,sd=diff(range(Design[,coord]))/500))
        if(!length(options$par))
          options$par<-NULL
        tempMax<-getMax(x=Design[i,coord],f = f,fprime = fprime,coord,d,options=options)
        results$max[i,coord]<-tempMax$val
        allMaxPoints[(coord-1)*nPtsPerDim+i,-coord]<-tempMax$aux$par
        allMaxPoints[(coord-1)*nPtsPerDim+i,coord]<-Design[i,coord]

        options$par<-allMinPoints[getClosePoints(x = Design[i,coord],allPts = allMinPoints,whichDim = coord),-coord]
        # 'regularization'
        options$par<-options$par*(1+rnorm(n=1,sd=diff(range(Design[,coord]))/500))
        if(!length(options$par))
          options$par<-NULL
        tempMin<-getMin(Design[i,coord],f = f,fprime = fprime,coord= coord,d = d,options=options)
        results$min[i,coord]<-tempMin$val
        allMinPoints[(coord-1)*nPtsPerDim+i,-coord]<-tempMin$aux$par
        allMinPoints[(coord-1)*nPtsPerDim+i,coord]<-Design[i,coord]

      }

      if(!is.null(fprime)){
        gr<-fprime(allMaxPoints[(coord-1)*nPtsPerDim+i,])
        gradOpt[i,coord,2]<-gr[coord]
        normGradRes[i,coord,2]<-norm(gr[-coord],"2")

        gr<-fprime(allMinPoints[(coord-1)*nPtsPerDim+i,])
        gradOpt[i,coord,1]<-gr[coord]
        normGradRes[i,coord,1]<-norm(gr[-coord],"2")
      }
    }

    if(!is.null(options$plts))
      if(options$plts){
        plot(Design[,coord],results$min[,coord],ylim=c(min(results$min[,coord]),max(results$max[,coord])),type='l',main=paste("Coordinate",coord))
        lines(Design[,coord],results$max[,coord])

 #       boxplot(tempMCvals[[coord]])
 #       lines(apply(tempMCvals[[coord]],2,max),col=2)
 #       lines(apply(tempMCvals[[coord]],2,min),col=2)
      }
    if(!is.null(options$verb))
      if(options$verb)
        cat("###############################\n")
  }

  toReturn<-list(res=results)
  if(is.null(options$heavyReturn))
    options$heavyReturn=FALSE
  if(options$heavyReturn){
    toReturn$minima<-allMinPoints
    toReturn$maxima<-allMaxPoints
    toReturn$gradRes<-normGradRes
    toReturn$gradOpt<-gradOpt
  }
  return(toReturn)
}



# approxMaxMin
#' @author Dario Azzimonti
#' @name approxMaxMin
#' @title Approximate coordinate profile functions
#' @description Evaluate profile extrema over other variables with approximations at few values
#' @param f the function to be evaluated
#' @param fprime derivative of the function
#' @param d dimension of the input domain
#' @param opts a list containing the options for this function and the subfunctions getMax, getMin or getMaxMinMC, see documentation of getMax, getMin, getMaxMinMC for details. The options only for approxMaxMin are
#' \itemize{
#' \item{\code{limits:}}{an optional list with the upper and lower limits of each dimension, if NULL then for each dimension limits are 0,1}
#' \item{\code{smoother:}}{Select which smoother to use:a string that selects which smoother to use: \itemize{
#'       \item{\code{"1order"}}: first order interpolation with gradient
#'       \item{\code{"splineSmooth"}}: smoothing spline with default degrees of freedom (DEFAULT OPTION)
#'       \item{\code{"quantSpline"}}: profile inf and profile sup approximated with quantile spline regression at levels 0.1 and 0.9 respectively
#' }}
#' \item{\code{heavyReturn:}}{If TRUE returns also all minimizers, default is FALSE.}
#' \item{\code{initDesign:}}{The design of few points where the expensive sup is evaluated.}
#' \item{\code{fullDesignSize:}}{The full design where the function is approximated.}
#' \item{\code{multistart:}}{number of multistarts for optim procedure.}
#' \item{\code{MonteCarlo:}}{if TRUE, computes sup with Monte Carlo procedure.}
#' \item{\code{numMCsamples:}}{number of MC samples for the sup.}
#' \item{\code{plts:}}{If TRUE, plots the max/min functions at each coordinate, default is FALSE.}
#' \item{\code{verb:}}{If TRUE, outputs intermediate results, default is FALSE.}
#' }
#' @return a list of two data frames (min, max) of the evaluations of \eqn{f_sup(x_i) = sup_{x_j \neq i} f(x_1,\dots,x_d) } and \eqn{f_inf(x_i) = inf_{x_j \neq i} f(x_1,\dots,x_d) }
#' for each i at the design Design. By default Design is a 100 equally spaced points for each dimension. It can be changed by defining it in options$Design
#' @examples
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # Compute the coordinate profile extrema with full optimization on 2d example
#'
#' # Define the function
#' g=function(x){
#'   return(-branin(x))
#' }
#' # Define the gradient
#' gprime = function(x){
#'   x1 = x[1]*15-5
#'   x2 = x[2]*15
#'   f1prime = (15*25)/(4*pi^4)*x1^3 - (15*75)/(2*pi^3)*x1^2 +
#'   (80*15)/(pi^2)*x1 - (5*15)/(pi^2)*x2*x1 +
#'   10*15/pi*x2 - 60*15/pi-10*15* (1 - 1/(8*pi))*sin(x1)
#'   f2prime = 2*15*(x2-5/(4*pi^2)*x1^2 +5/pi*x1-6)
#'   return(matrix(c(-f1prime,-f2prime),nrow=1))
#' }
#'
#' # generic approximation options
#' init_des<-lhs::maximinLHS(15,2)
#' options_approx<- list(multistart=4,heavyReturn=TRUE,initDesign=init_des,fullDesignSize=100)
#'
#' # 1order approximation
#' options_approx$smoother<-"1order"
#' coordProf_approx_1order<-approxMaxMin(f = g,fprime = gprime,d=2,opts = options_approx)
#'
#' # quantile regression
#' options_approx$smoother<-"quantSpline"
#' coordProf_approx_quantReg<-approxMaxMin(f = g,fprime = gprime,d=2,opts = options_approx)
#'
#'
#' \donttest{
#' # Consider threshold=-10
#' threshold<- -10
#' # obtain the points where the profiles take the threshold value
#' pp_change<-getChangePoints(threshold = threshold,allRes = coordProf_approx_quantReg)
#' # evaluate g at a grid and plot the image
#' x<-seq(0,1,,100)
#' grid<-expand.grid(x,x)
#' g_evals<- apply(X = grid,MARGIN = 1,FUN = g)
#' image(x = x,y = x,z = matrix(g_evals,nrow = 100),col = grey.colors(20))
#' contour(x=x,y=x,z=matrix(g_evals,nrow = 100), add=TRUE, nlevels = 20)
#' contour(x=x,y=x,z=matrix(g_evals,nrow = 100), add=TRUE, levels = threshold,col=2)
#' abline(h = pp_change$neverEx$`-10`[[2]],col="darkgreen",lwd=2)
#' abline(v = pp_change$neverEx$`-10`[[1]],col="darkgreen",lwd=2)
#' # Plot the coordinate profiles and a threshold
#' plotMaxMin(allRes = coordProf_approx_1order,threshold = threshold,changes = TRUE)
#' plotMaxMin(allRes = coordProf_approx_quantReg,threshold = threshold,changes = TRUE)
#' }
#' @export
approxMaxMin=function(f,fprime=NULL,d,opts=NULL){
  if(is.null(opts$limits)){
    limits<-list(lower=rep(0,d),upper=rep(1,d))
  }
#  nPtsPerDim<-rep(100,d)
  if(is.null(opts$initDesign)){
    opts$initDesign<-maximinLHS(ceiling(sqrt(d)*10),d)
  }
  # It is better to include points on the border
#  if(!is.null(opts$smoother) && opts$smoother=="quantSpline"){
    for(coord in seq(d)){
      opts$initDesign[,coord]<-c(limits$lower[coord],sort(opts$initDesign[,coord])[c(-1,-nrow(opts$initDesign))],limits$upper[coord])
    }
#  }
  if(is.null(opts$fullDesignSize)){
    opts$fullDesignSize<-80*d
  }

  aa<-getAllMaxMin(f,fprime,d,options=list(Design=opts$initDesign,multistart=opts$multistart,
                                           heavyReturn=T,MonteCarlo=opts$MonteCarlo,numMCsamples=opts$numMCsamples,verb=opts$verb))

  # initialize the max/min values
  sK_max<-data.frame(matrix(NA,nrow=opts$fullDesignSize,ncol=d))
  sK_min<-data.frame(matrix(NA,nrow=opts$fullDesignSize,ncol=d))

  if(!is.null(fprime)){
    aa$res$grad_min<-matrix(NA,ncol=ncol(aa$res$min),nrow=nrow(aa$res$min))
    colnames(aa$res$grad_min)<-colnames(aa$res$min)

    aa$res$grad_max<-matrix(NA,ncol=ncol(aa$res$max),nrow=nrow(aa$res$max))
    colnames(aa$res$grad_max)<-colnames(aa$res$max)
  }

  for(coord in seq(d)){
    # the gradient of the sup/inf function is saved
    nn<-length(aa$res$max[,coord])
    for(i in seq(nn)){
      if(!is.null(fprime)){
        yy<-aa$minima[(coord-1)*nn+i,]
        aa$res$grad_min[i,coord]<-fprime(yy)[,coord]
        yy<-aa$maxima[(coord-1)*nn+i,]
        aa$res$grad_max[i,coord]<-fprime(yy)[,coord]
      }

    }
    newPoints<-matrix(seq(limits$lower[coord],limits$upper[coord],,opts$fullDesignSize),ncol=1)

    if(!is.null(opts$smoother) && opts$smoother=="1order"){
      sK_max[,coord]<- kGradSmooth(newPoints=newPoints,profPoints = aa$maxima[((coord-1)*nn+1):(coord*nn),coord],
                                   profEvals=aa$res$max[,coord], profGradient=aa$res$grad_max[,coord])

      sK_min[,coord]<-kGradSmooth(newPoints=newPoints,profPoints = aa$minima[((coord-1)*nn+1):(coord*nn),coord],
                                  profEvals=aa$res$min[,coord], profGradient=aa$res$grad_min[,coord])
    }else if(!is.null(opts$smoother) && opts$smoother=="quantSpline"){
      dd_fit<-data.frame(x=aa$maxima[((coord-1)*nn+1):(coord*nn),coord],y=aa$res$max[,coord])
      fit.max <- rq(y ~ bs(x, df=max(nn%/%1.5,3)),data=dd_fit,tau=0.9)
      dd_fit<-data.frame(x=aa$minima[((coord-1)*nn+1):(coord*nn),coord],y=aa$res$min[,coord])
      fit.min <- rq(y ~ bs(x, df=max(nn%/%1.5,3)), tau=0.1,data=dd_fit)
      newPoints<-data.frame(x=newPoints)

      sK_max[,coord]<- predict.rq(object = fit.max,newdata = newPoints)
      sK_min[,coord]<- predict.rq(object = fit.min,newdata = newPoints)
    }else {
      sk_maxSmooth<-smooth.spline(x=aa$maxima[((coord-1)*nn+1):(coord*nn),coord],y=aa$res$max[,coord])
      sk_minSmooth<-smooth.spline(x=aa$minima[((coord-1)*nn+1):(coord*nn),coord],y=aa$res$min[,coord])

      sK_max[,coord]<-predict(sk_maxSmooth,newPoints)$y
      sK_min[,coord]<-predict(sk_minSmooth,newPoints)$y
    }




    ## Debug plots
    if(!is.null(opts$debug) && opts$debug){
  #    par(mfrow=c(1,2))
      # Plot  max/min function
      ylimTemp<-range(c(sK_min[,coord],aa$res$min[,coord],sK_max[,coord],aa$res$max[,coord]))
      plot(aa$minima[((coord-1)*nn+1):(coord*nn),coord],aa$res$min[,coord],main=paste("Coordinate",coord,", Min"),
           xlab="x",ylab="f",ylim=ylimTemp)
      lines(newPoints,sK_min[,coord],type='l',col=3)
      points(aa$maxima[((coord-1)*nn+1):(coord*nn),coord],aa$res$max[,coord])
      lines(newPoints,sK_max[,coord],type='l',col=3)

      if(!is.null(opts$smooth) && opts$smooth){
        points(sk_minSmooth$x,sk_minSmooth$y,col=2)
        points(sk_maxSmooth$x,sk_maxSmooth$y,col=2)
      }


      # Plot min function
#      ylimTemp<-range(c(sK_min[,coord],aa$res$min[,coord]))
#      plot(aa$minima[((coord-1)*nn+1):(coord*nn),coord],aa$res$min[,coord],main=paste("Coordinate",coord,", Min"),
#           xlab="x",ylab="f")
#      lines(newPoints,sK_min[,coord],type='l',col=3)
#      points(sk_minSmooth$x,sk_minSmooth$y,col=2)

      # This code uses gg_min, which needs to be initialized before the function
#      ylimTemp<-c(min(c(gg_min[,coord],aa$res$min[,coord])),max(c(gg_min[,coord],aa$res$min[,coord])))
#      plot(dd,gg_min[,coord],ylim=ylimTemp,type='l',main=paste("Coordinate",coord,", Min"),
#           xlab="x",ylab="f")
#      points(aa$minima[((coord-1)*nn+1):(coord*nn),coord],aa$res$min[,coord])
#      points(aa$minima[((coord-1)*nn+1):(coord*nn),coord],aa$res$grad_min[,coord],pch=17)
#      lines(dd,allRes_km$res$min[,coord])
#      abline(h=0,col=3)
#      abline(h=-10,col=2)
#      lines(newPoints,sK_min[,coord],col=4)
      #  points(allMinPoints[(coord-1)*nPtsPerDim[coord]+1:i,coord],results$min[1:i,coord])

      # Plot the max function
#      ylimTemp<-range(c(sK_max[,coord],aa$res$max[,coord]))
#      plot(aa$maxima[((coord-1)*nn+1):(coord*nn),coord],aa$res$max[,coord],main=paste("Coordinate",coord,", Max"),
#           xlab="x",ylab="f")
#      lines(newPoints,sK_max[,coord],type='l',col=3)
    #  lines(pp_sk_max_smooth$x,pp_sk_max_smooth$y,type='l',col=2)
#      points(sk_maxSmooth$x,sk_maxSmooth$y,col=2)

      # This code uses gg_min, which needs to be initialized before the function
#      ylimTemp<-c(min(c(gg_max[,coord],aa$res$max[,coord])),max(c(gg_max[,coord],aa$res$max[,coord])))
#      plot(dd,gg_max[,coord],ylim=ylimTemp,type='l',main=paste("Coordinate",coord,", Max"),
#           xlab="x",ylab="f")
#      points(aa$maxima[((coord-1)*nn+1):(coord*nn),coord],aa$res$max[,coord])
#      points(aa$maxima[((coord-1)*nn+1):(coord*nn),coord],aa$res$grad_max[,coord],pch=17)
#      lines(dd,allRes_km$res$max[,coord])
#      abline(h=0,col=3)
#      abline(h=-10,col=2)
#      lines(newPoints,sK_max[,coord],col=4)
      #  points(allMaxPoints[(coord-1)*nPtsPerDim[coord]+1:i,coord],results$max[1:i,coord])
    }
  }

  res=list(min=sK_min,max=sK_max)
  toReturn<-list(res=res)

  if(!is.null(opts$heavyReturn) && opts$heavyReturn){
    aa$design<-opts$initDesign
    toReturn$profPoints<-aa
  }
  return(toReturn)

}


#' @author Dario Azzimonti
#' @title First order approximation
#' @name kGradSmooth
#' @description Compute first order approximation of function from evaluations and gradient
#' @param newPoints vector of points where to approximate the function
#' @param profPoints locations where the function was evaluated
#' @param profEvals value of the evaluation at profPoints
#' @param profGradient value of the gradient at profPoints
#' @return approximated values of the function at newPoints
kGradSmooth<-function(newPoints,profPoints,profEvals,profGradient){
  newPoints<-matrix(newPoints,ncol=1)
  clPts<-apply(newPoints,MARGIN = 1,FUN = function(x){which.min(abs(x-profPoints))})
  ww<-rep(NA,nrow(newPoints))
  for(i in seq(nrow(newPoints))){
    ww[i]<-profEvals[clPts[i]]+(newPoints[i,]-profPoints[clPts[i]])*profGradient[clPts[i]]
  }
  return(ww)
}


#' @author Dario Azzimonti
#' @title Find close points
#' @name getClosePoints
#' @description Obtain points close in one specific dimension
#' @param x one dimensional point
#' @param allPts dataframe containing a list of d dimensional points
#' @param whichDim integer defining the dimension of x
#' @return the index in allPts (row number) of the closest point in allPts to x along the whichDim dimension
getClosePoints<-function(x,allPts,whichDim){
  rr<-which.min(abs(x-allPts[,whichDim]))
  if(!length(rr))
    rr<-NULL
  return(rr)
}
