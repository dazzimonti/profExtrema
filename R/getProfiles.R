# getProfileSup function
#' @author Dario Azzimonti
#' @name getProfileSup
#' @title Generic profile sup function computation
#' @description Compute profile sup function for an arbitrary matrix \code{Psi} with \link[nloptr]{nloptr}.
#' @param eta one dimensional point where the function is to be evaluated
#' @param Psi projection matrix
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format)
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to optim:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d-1)}
#' \item{\code{lower:} }{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{max_{x_1,\dots, x_{coord-1},x_{coord+1}, \dots, x_d} f(x_1,\dots,x_d)}
#' @export
# requires library(nloptr)
getProfileSup = function(eta,Psi,f,fprime,d,options=NULL){

  lower=rep(0,d)
  upper=rep(1,d)
  trace=0

  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace

  # Verify if the feasible set is not empty
  Hlarge<-makeH(a1 = rbind(diag(-1,nrow = d,ncol=d),diag(1,nrow = d,ncol=d)),
                b1 = rbind(matrix(lower-1e-8,nrow = d),matrix(upper-1e-8,nrow = d)),
                a2 = Psi,b2=matrix(eta))
  solLP<-lpcdd(hrep = Hlarge,objgrd = rep(0,d))

  if(solLP$solution.type=="Optimal"){
    # Find any solution of the linear system
    startingPoint<-matrix(solLP$primal.solution,nrow=d)
  }else{
    warning("Cannot proceed with these constraints. A simple LP optimization problem does not have an optimal solution")
    optRes<-list(value=NA,solLP=solLP,message="The feasible set is either empty (solLP$solution.type==Inconsistent) or the problem is dual inconsistent. See solLP for the results of a LP problem with those contraints and objective function 0.")
    return(list(val=optRes$value,aux=optRes))
  }


  # nloptr does minimizations so we need to change the sign
  ff<-function(x){return(-f(x))}
  ffprime<-function(x){return(-fprime(x))}




  if(!is.null(options$par)){
    startingPoint<-options$par
    if(Psi%*%startingPoint< eta)
      warning("Provided starting point does not respect constraints")
  }
#  }else{
#    startingPoint<-runif(min = lower,max = upper,d)
#  }

  eval_g0 <-function(x){
    return(Psi%*%x -eta)
  }

  eval_jac_g0 <- function(x){
    return(Psi)
  }

  if(is.null(options$opts)){
    local_opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                        "xtol_rel"  = 1.0e-6 )
    opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                  "xtol_rel"  = 1.0e-6,
                  "maxeval"   = 1000,
                  "local_opts" = local_opts )
  }else{
    opts <- options$opts
  }
  if(!is.null(options$multistart)){
  #  cat("not implemented yet!")
    optRess<-NULL
    optRes<-list(objective=-Inf,aux=NA)
    for(i in seq(options$multistart)){
      if(i>1)
        startingPoint<-runif(d,min = lower,max=upper)

      optRess<-nloptr(x0=startingPoint, eval_f=ff,
                      eval_grad_f=ffprime, eval_g_eq = eval_g0,
                      eval_jac_g_eq = eval_jac_g0,
                      opts=opts,lb=lower,ub=upper)
      optRess$objective= -optRess$objective
      #      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
      if(optRess$objective>=optRes$objective){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    optRes<-nloptr(x0=startingPoint, eval_f=ff,
                   eval_grad_f=ffprime, eval_g_eq = eval_g0,
                   eval_jac_g_eq = eval_jac_g0,
                   opts=opts,lb=lower,ub=upper)
    optRes$objective=-optRes$objective
  }

  return(list(val=optRes$objective,aux=optRes))
}


# getProfileInf function
#' @author Dario Azzimonti
#' @name getProfileInf
#' @title Generic profile inf function computation
#' @description Compute profile inf function for an arbitrary matrix \code{Psi} with \link[nloptr]{nloptr}.
#' @param eta one dimensional point where the function is to be evaluated
#' @param Psi projection matrix
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format)
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to optim:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d)}
#' \item{\code{lower:} }{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{max_{x_1,\dots, x_{coord-1},x_{coord+1}, \dots, x_d} f(x_1,\dots,x_d)}
#' @export
# requires library(nloptr)
getProfileInf = function(eta,Psi,f,fprime,d,options=NULL){

  lower=rep(0,d)
  upper=rep(1,d)
  trace=0

  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace


  if(!is.null(options$par)){
    startingPoint<-options$par
    if(Psi%*%startingPoint< eta)
      warning("Provided starting point does not respect constraints")
  }else{
    startingPoint<-runif(min = lower,max = upper,d)
  }

  eval_g0 <-function(x){
    return(Psi%*%x -eta)
  }

  eval_jac_g0 <- function(x){
    return(Psi)
  }

  if(is.null(options$opts)){
    local_opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                        "xtol_rel"  = 1.0e-6 )
    opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                  "xtol_rel"  = 1.0e-6,
                  "maxeval"   = 1000,
                  "local_opts" = local_opts )
  }else{
    opts <- options$opts
  }

  if(!is.null(options$multistart)){
    optRess<-NULL
    optRes<-list(objective=Inf,aux=NULL)
    for(i in seq(options$multistart)){
      startingPoint<-runif(d,min = lower,max=upper)
      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      optRess<-nloptr(x0=startingPoint, eval_f=f,
                     eval_grad_f=fprime, eval_g_eq = eval_g0,
                     eval_jac_g_eq = eval_jac_g0,
                     opts=opts,lb=lower,ub=upper)
#      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
      if(optRess$objective<optRes$objective){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    optRes<-nloptr(x0=startingPoint, eval_f=f,
                   eval_grad_f=fprime, eval_g_eq = eval_g0,
                   eval_jac_g_eq = eval_jac_g0,
                   opts=opts,lb=lower,ub=upper)
  }

  return(list(val=optRes$objective,aux=optRes))
}


#' @author Dario Azzimonti
#' @name getProfileCoord
#' @title Coordinate profile extrema with nloptr
#' @description Compute coordinate profile extrema with \link[nloptr]{nloptr}
#' @param f the function to be evaluated
#' @param fprime derivative of the function
#' @param d dimension of the input domain
#' @param options a list containing the options for this function and the subfunctions \link{getMaxMinMC},
#' see documentation of \link{getProfileSup}, \link{getProfileInf} for details. The options only for getProfileCoord are
#' \itemize{
#' \item{\code{Design:}}{an optional design matrix with the discretization of each dimension, if NULL then for each dimension Design[,coord] = seq(0,1,length.out=100)}
#' \item{\code{heavyReturn:}}{If TRUE returns also all minimizers, default is FALSE.}
#' \item{\code{plts:}}{If TRUE, plots the max/min functions at each coordinate, default is FALSE.}
#' \item{\code{verb:}}{If TRUE, outputs intermediate results, default is FALSE.}
#' \item{\code{MonteCarlo:}}{If TRUE, use the MC optimizer otherwise use standard optim, default is FALSE.}
#' }
#' @return a list of two data frames (min, max) of the evaluations of \eqn{f_sup(x_i) = sup_{x_j \neq i} f(x_1,\dots,x_d) } and \eqn{f_inf(x_i) = inf_{x_j \neq i} f(x_1,\dots,x_d) }
#' for each i at the design \code{Design}. By default \code{Design} is a 100 equally spaced points for each dimension. It can be changed by defining it in \code{options$Design}
#' @export
# REUSE THE MC SAMPLES? you need to save them find a way to find the closest point to the evaluated one
# and then check if the values are consistent?
getProfileCoord<-function(f,fprime=NULL,d,options=NULL){

  if(is.null(options$lower))
    options$lower=rep(0,d)

  if(is.null(options$upper))
    options$upper=rep(1,d)

  if(is.null(options$Design)){
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
    Psi<-matrix(rep(0,d),ncol=d)
    Psi[,coord]<-1
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
        options$par<-allMaxPoints[getClosePoints(x = Design[i,coord],allPts = allMaxPoints,whichDim = coord),]
        options$par<-options$par*(1-runif(n=1,min=0,max=diff(range(Design[,coord]))/500))#(1+rnorm(n=1,sd=diff(range(Design[,coord]))/500))
        if(!length(options$par)){
          options$par<-NULL
        }else{
          options$par[coord]<-Design[i,coord]
        }
        tempMax<-getProfileSup(eta = Design[i,coord],Psi = Psi,f = f,fprime = fprime,d = d,options=options)
        results$max[i,coord]<-tempMax$val
        allMaxPoints[(coord-1)*nPtsPerDim+i,]<-tempMax$aux$solution
#        allMaxPoints[(coord-1)*nPtsPerDim+i,coord]<-Design[i,coord]

        options$par<-allMinPoints[getClosePoints(x = Design[i,coord],allPts = allMinPoints,whichDim = coord),]
        # 'regularization'
        options$par<-options$par*(1-runif(n=1,min=0,max=diff(range(Design[,coord]))/500))
        if(!length(options$par)){
          options$par<-NULL
        }else{
          options$par[coord]<-Design[i,coord]
        }
        tempMin<-getProfileInf(eta=Design[i,coord],Psi=Psi,f = f,fprime = fprime,d = d,options=options)
        results$min[i,coord]<-tempMin$val
        allMinPoints[(coord-1)*nPtsPerDim+i,]<-tempMin$aux$solution
  #      allMinPoints[(coord-1)*nPtsPerDim+i,coord]<-Design[i,coord]

      }

      gr<-fprime(allMaxPoints[(coord-1)*nPtsPerDim+i,])
      gradOpt[i,coord,2]<-gr[coord]
      normGradRes[i,coord,2]<-norm(gr[-coord],"2")

      gr<-fprime(allMinPoints[(coord-1)*nPtsPerDim+i,])
      gradOpt[i,coord,1]<-gr[coord]
      normGradRes[i,coord,1]<-norm(gr[-coord],"2")
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

# approxCoordProf
#' @author Dario Azzimonti
#' @name approxCoordProf
#' @title Approximate coordinate profiles with nloptr
#' @description Evaluate profile extrema over other variables with approximations at few values
#' @param f the function to be evaluated
#' @param fprime derivative of the function
#' @param optimFun function to do full optimization: at the moment either getAllMaxMin (default choice) or getProfileCoord
#' @param threshold threshold of interest
#' @param d dimension of the input domain
#' @param opts a list containing the options for this function and the subfunctions getMax, getMin or getMaxMinMC, see documentation of getMax, getMin, getMaxMinMC for details. The options only for approxMaxMin are
#' \itemize{
#' \item{\code{limits:}}{an optional list with the upper and lower limits of each dimension, if NULL then for each dimension limits are 0,1}
#' \item{\code{smoother:}}{a string that selects which smoother to use: \itemize{
#'       \item{\code{"1order"}}: first order interpolation with gradient
#'       \item{\code{"splineSmooth"}}: smoothing spline with default df (DEFAULT OPTION)
#'       \item{\code{"quantSpline"}}: quantile spline regression at level 0.95
#' }}
#' \item{\code{heavyReturn:}}{If TRUE returns also all minimizers, default is FALSE.}
#' \item{\code{initDesign:}}{The design of few points where the expensive sup is evaluated.}
#' \item{\code{fullDesignSize:}}{The full design where the function is approximated.}
#' \item{\code{preWork:}}{If TRUE, runs a pre-processing procedure to adjust the bad points.}
#' \item{\code{multistart:}}{number of multistarts for optim procedure.}
#' \item{\code{MonteCarlo:}}{if TRUE, computes sup with Monte Carlo procedure.}
#' \item{\code{numMCsamples:}}{number of MC samples for the sup.}
#' \item{\code{plts:}}{If TRUE, plots the max/min functions at each coordinate, default is FALSE.}
#' \item{\code{verb:}}{If TRUE, outputs intermediate results, default is FALSE.}
#' }
#' @return a list of two data frames (min, max) of the evaluations of \eqn{f_sup(x_i) = sup_{x_j \neq i} f(x_1,\dots,x_d) } and \eqn{f_inf(x_i) = inf_{x_j \neq i} f(x_1,\dots,x_d) }
#' for each i at the design Design. By default Design is a 100 equally spaced points for each dimension. It can be changed by defining it in options$Design
approxProfCoord=function(f,fprime=NULL,optimFun=getAllMaxMin,threshold=0,d,opts=NULL){
  if(is.null(opts$limits)){
    limits<-list(lower=rep(0,d),upper=rep(1,d))
  }
  #  nPtsPerDim<-rep(100,d)
  if(is.null(opts$initDesign)){
    opts$initDesign<-maximinLHS(ceiling(sqrt(d)*10),d)
  }
  if(is.null(opts$fullDesignSize)){
    opts$fullDesignSize<-80*d
  }

  aa<-optimFun(f,fprime,d,options=list(Design=opts$initDesign,multistart=opts$multistart,
                                           heavyReturn=T,MonteCarlo=opts$MonteCarlo,numMCsamples=opts$numMCsamples,verb=opts$verb))

  # initialize the max/min values
  sK_max<-data.frame(matrix(NA,nrow=opts$fullDesignSize,ncol=d))
  sK_min<-data.frame(matrix(NA,nrow=opts$fullDesignSize,ncol=d))


  # prework to adjust the ``bad'' points
  # if(!is.null(opts$preWork) && opts$preWork){
  #   # For each point in minima/maxima close to threshold let's ``resample''
  #   closeT<-threshold*(1-sign(threshold)*0.05)
  #   numOfresamples<-max(c(apply(aa$res$max>closeT,2,sum),apply(aa$res$min>closeT,2,sum)))
  #
  #   addResults<-matrix(NA,ncol=d,nrow=numOfresamples)
  #   colnames(addResults)<-colnames(aa$res$max)
  #   for (coord in seq(d)){
  #     reS_max<-c(opts$initDesign[which(aa$res$max[,coord]>closeT),coord], runif(numOfresamples-length(which(aa$res$max[,coord]>closeT))))
  #     reS_min<-c(opts$initDesign[which(aa$res$min[,coord]>closeT),coord], runif(numOfresamples-length(which(aa$res$min[,coord]>closeT))))
  #
  #     newS_max=reS_max*(1+rnorm(numOfresamples,sd=(limits$upper-limits$lower)[coord]/50))
  #     newS_min=reS_min*(1+rnorm(numOfresamples,sd=(limits$upper-limits$lower)[coord]/50))
  #
  #     #    plot(reS_min,newS_min)
  #     #    plot(reS_max,newS_max)
  #     # Make space for the new results
  #     # make space in the max/min results
  #     if(coord==1){
  #       aa$res$max<-rbind(aa$res$max,addResults)
  #       aa$res$min<-rbind(aa$res$min,addResults)
  #     }
  #     nn<-nrow(opts$initDesign)
  #
  #     # make space in the maxima/minima results
  #     aa$maxima<-rbind(aa$maxima[(1):(nn*coord+numOfresamples*(coord-1)),],addResults,
  #                      aa$maxima[-((1):(nn*coord+numOfresamples*(coord-1))),])
  #     aa$minima<-rbind(aa$minima[(1):(nn*coord+numOfresamples*(coord-1)),],addResults,
  #                      aa$minima[-((1):(nn*coord+numOfresamples*(coord-1))),])
  #
  #     for(i in seq(numOfresamples)){
  #
  #       tempMax<-getMax(newS_max[i],f = f,fprime = fprime,coord,d)
  #       aa$res$max[i+nrow(opts$initDesign),coord]<-tempMax$val
  #       aa$maxima[(coord)*nrow(opts$initDesign)+i+(coord-1)*numOfresamples,-coord]<-tempMax$aux$par
  #       aa$maxima[(coord)*nrow(opts$initDesign)+i+(coord-1)*numOfresamples,coord]<-newS_max[i]
  #
  #       par(mfrow=c(1,2))
  #       ylimTemp<-c(min(c(allRes_km$res$max[,coord],aa$res$max[,coord],threshold),na.rm = T),max(c(allRes_km$res$max[,coord],aa$res$max[,coord],threshold),na.rm = T))
  #       plot(dd,allRes_km$res$max[,coord],ylim=ylimTemp,type='l',main=paste("Coordinate",coord,", Max, pp",i),
  #            xlab="x",ylab="f")
  #       points(aa$maxima[((coord-1)*(numOfresamples+nrow(opts$initDesign))+1):((coord-1)*(numOfresamples)+coord*nrow(opts$initDesign)+i),coord],aa$res$max[(1):(nrow(opts$initDesign)+i),coord])
  #       points(aa$maxima[(coord)*nrow(opts$initDesign)+i+(coord-1)*numOfresamples,coord],aa$res$max[i+nrow(opts$initDesign),coord],pch=17)
  #       #  abline(h=0,col=3)
  #       abline(h=-10,col=2)
  #
  #       if(i<length(which(aa$res$max[,coord]>closeT))){
  #         indRecompute<-which(aa$res$max[,coord]>closeT)[i]
  #         tempMax<-getMax(reS_max[i],f = f,fprime = fprime,coord,d,options=list(par=getClosePoints(aa$maxima[indRecompute,-coord],aa$maxima,whichDim = coord)))
  #         aa$res$max[indRecompute,coord]<-max(tempMax$val,aa$res$max[indRecompute,coord])
  #         aa$maxima[indRecompute+(coord-1)*(nrow(opts$initDesign)+numOfresamples),-coord]<-tempMax$aux$par
  #         aa$maxima[indRecompute+(coord-1)*(nrow(opts$initDesign)+numOfresamples),coord]<-reS_max[i]
  #
  #         points(reS_max[i],aa$res$max[indRecompute,coord],pch=16,col=2)
  #       }
  #
  #
  #
  #       tempMin<-getMin(newS_min[i],f = f,fprime = fprime,coord,d)
  #       aa$res$min[i+nrow(opts$initDesign),coord]<-tempMin$val
  #       aa$minima[(coord)*nrow(opts$initDesign)+i+(coord-1)*numOfresamples,-coord]<-tempMin$aux$par
  #       aa$minima[(coord)*nrow(opts$initDesign)+i+(coord-1)*numOfresamples,coord]<-newS_min[i]
  #
  #       ylimTemp<-c(min(c(allRes_km$res$min[,coord],aa$res$min[,coord],threshold),na.rm = T),max(c(allRes_km$res$min[,coord],aa$res$min[,coord],threshold),na.rm = T))
  #       plot(dd,allRes_km$res$min[,coord],ylim=ylimTemp,type='l',main=paste("Coordinate",coord,", Min, pp",i),
  #            xlab="x",ylab="f")
  #       points(aa$minima[((coord-1)*(numOfresamples+nrow(opts$initDesign))+1):((coord-1)*(numOfresamples)+coord*nrow(opts$initDesign)+i),coord],aa$res$min[(1):(nrow(opts$initDesign)+i),coord])
  #       points(aa$minima[(coord)*nrow(opts$initDesign)+i+(coord-1)*numOfresamples,coord],aa$res$min[i+nrow(opts$initDesign),coord],pch=17)
  #       #  abline(h=0,col=3)
  #       abline(h=-10,col=2)
  #
  #       if(i<length(which(aa$res$min[,coord]>closeT))){
  #         indRecompute<-which(aa$res$min[,coord]>closeT)[i]
  #         tempMin<-getMin(reS_min[i],f = f,fprime = fprime,coord,d,options=list(par=getClosePoints(aa$maxima[indRecompute,-coord],aa$maxima,whichDim = coord)))
  #         aa$res$min[indRecompute,coord]<-min(tempMin$val,aa$res$min[indRecompute,coord])
  #         aa$minima[indRecompute+(coord-1)*(nrow(opts$initDesign)+numOfresamples),-coord]<-tempMin$aux$par
  #         aa$minima[indRecompute+(coord-1)*(nrow(opts$initDesign)+numOfresamples),coord]<-reS_min[i]
  #
  #         points(reS_min[i],aa$res$min[indRecompute,coord],pch=16,col=2)
  #       }
  #       invisible(readline(prompt="Press [enter] to continue"))
  #
  #     }
  #   }
  # }

  aa$res$grad_min<-matrix(NA,ncol=ncol(aa$res$min),nrow=nrow(aa$res$min))
  colnames(aa$res$grad_min)<-colnames(aa$res$min)

  aa$res$grad_max<-matrix(NA,ncol=ncol(aa$res$max),nrow=nrow(aa$res$max))
  colnames(aa$res$grad_max)<-colnames(aa$res$max)

  for(coord in seq(d)){
    # the gradient of the sup/inf function is saved
    nn<-length(aa$res$max[,coord])
    for(i in seq(nn)){
      yy<-aa$minima[(coord-1)*nn+i,]
      aa$res$grad_min[i,coord]<-matrix(fprime(yy),ncol=d)[,coord]

      yy<-aa$maxima[(coord-1)*nn+i,]
      aa$res$grad_max[i,coord]<-matrix(fprime(yy),ncol=d)[,coord]
    }
    newPoints<-matrix(seq(limits$lower[coord],limits$upper[coord],,opts$fullDesignSize),ncol=1)

    if(!is.null(opts$smoother) && opts$smoother=="1order"){
      sK_max[,coord]<- kGradSmooth(newPoints=newPoints,profPoints = aa$maxima[((coord-1)*nn+1):(coord*nn),coord],
                                   profEvals=aa$res$max[,coord], profGradient=aa$res$grad_max[,coord])

      sK_min[,coord]<-kGradSmooth(newPoints=newPoints,profPoints = aa$minima[((coord-1)*nn+1):(coord*nn),coord],
                                  profEvals=aa$res$min[,coord], profGradient=aa$res$grad_min[,coord])
    }else if(opts$smoother=="quantSpline"){
      fit.max <- rq(aa$res$max[,coord] ~ bs(aa$maxima[((coord-1)*nn+1):(coord*nn),coord], df=15), tau=0.95)
      fit.min <- rq(aa$res$min[,coord] ~ bs(aa$minima[((coord-1)*nn+1):(coord*nn),coord], df=15), tau=0.05)
      newPoints<-data.frame(newPoints)

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
