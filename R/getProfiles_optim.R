# getProfileSup function with smarter algorithm!
#' @author Dario Azzimonti
#' @name solutionLinConstraint
#' @title Obtain generic solution to linear constraints
#' @description Computes an arbitrary solution \eqn{x^*} that satisfies \code{Phi x = eta}
#' @param eta \eqn{p} dimensional point where the function is to be evaluated
#' @param Phi \eqn{p x d} matrix describing the projection matrix
#' @param d dimension of the input
#' @param p dimension of the projected space
#' @param par a vector of dimension d-p that can be selected by the user to guide the solution. When NULL (suggested) it is chosen by the function.
#' @return a vector of dimensions \code{d x 1} containing the solution \eqn{x^*}.
# requires library(MASS) and lpSolve
solutionLinConstraint<- function(eta,Phi,d,p,par=NULL){
  kk<-0
  xiSol<-NULL
  allCombs<-combn(1:d,p)

  while(is.null(xiSol)){
    kk<-kk+1
    for(dd in seq(choose(d,p))){
      if(is.null(par)){
        if(kk==1){
          arbRest<-rep(0,d-p)
        }else if(kk==2){
          arbRest<-rep(1,d-p)
        }else{
          arbRest<-runif(d-p)
        }
      }else{
        arbRest=par
      }
      solvable<-Phi[1:p,allCombs[,dd]]
      rest<-matrix(Phi[1:p,setdiff(1:d,allCombs[,dd])],ncol=d-p)
      if(kappa(solvable)<1e16){
        tempSol<-matrix(c(solve(solvable,eta-rest%*%arbRest),arbRest))
        reOrderedSol<- tempSol
        reOrderedSol[allCombs[,dd]]<-tempSol[1:p]
        reOrderedSol[setdiff(1:d,allCombs[,dd])]<-tempSol[(p+1):d]
        reOrderedSol<-round(reOrderedSol,digits = 12)
      }else{
        next
      }

      if(all(reOrderedSol>=rep(0,d))&all(reOrderedSol<=rep(1,d))){
        if(is.null(xiSol))
          xiSol<-reOrderedSol
      }
    }
  }
  return(xiSol)
}


# getProfileSup function with smarter algorithm!
#' @author Dario Azzimonti
#' @name getProfileSup_optim
#' @title Generic profile sup function computation with optim
#' @description Compute profile sup function for an arbitrary matrix \code{Phi} with the L-BFGS-B algorithm of \link[stats]{optim}.
#' @param eta \eqn{p} dimensional point where the function is to be evaluated
#' @param Phi projection matrix of dimensions \code{p x d}
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format, returning a \eqn{d} dimensional vector)
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to optim:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d-1)}
#' \item{\code{lower:} }{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{max_{x \in D_Psi} f(x)}
#' @seealso \link{getProfileInf_optim}, \link{getProfileSup}, \link{plotMaxMin}
#' @export
# requires library(MASS) and lpSolve
getProfileSup_optim = function(eta,Phi,f,fprime,d,options=NULL){

  p = dim(Phi)[1]

  # put to zero small values
  Phi[which(abs(Phi)<.Machine$double.neg.eps)]<-0

  if(!is.null(options$trace))
    trace=options$trace

  # Find any solution of the linear system
  xiSol<-solutionLinConstraint(eta=eta,Phi=Phi,d = d,p=p)

  nullSpacePhi<-Null(t(Phi))

  # find equations for boundary conditions on z
  lower<-rep(0,d)
  upper<-rep(1,d)

  Az <- rbind(-nullSpacePhi,nullSpacePhi)
  bz <- rbind(xiSol-lower,-xiSol+upper)
  #  vertices_poly<-matrix(enumerate.vertices(A=Az,b=bz),ncol=d-p)
  rcdd_Hrep<-makeH(a1 = Az,b1 = bz)
  vertices_poly<-matrix(scdd(rcdd_Hrep)$output[,-c(1,2)],ncol=d-p)

  # generate point inside
  #  library(lpSolve)



  ff<-function(z){
    x_cand<-xiSol+nullSpacePhi%*%z
    g_cand<- sqrt(.Machine$double.xmax)
    if(all(x_cand>=rep(0,d))&all(x_cand<=rep(1,d))){
      return(-f(x_cand))
    }else{
      return(g_cand)
    }
  }

  ffprime<-function(z){
    x_cand<-xiSol+nullSpacePhi%*%z
    if(all(Az%*%z-bz<=0)){
      return(-crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))
    }else{
      ww<-Az%*%z-bz
      dirs<-which(ww>0)
      ww<-matrix(ww[dirs]/sum(ww[dirs]))
      return(-Az[dirs,]%*%ww)
    }
    return(-crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))
  }

  if(!is.null(options$par)){
    #    startingPoint<-options$par
    #   if(Phi%*%startingPoint< eta)
    warning("Constraints check not implemented!")
  }else{
    #    startingPoint<-runif(min = lower,max = upper,d-p)
    startingPoint<-colMeans(vertices_poly)#lp(direction = "max",objective.in = runif(d-p),const.mat = Az,const.dir = "<=",const.rhs = bz)$solution
    #      matrix(c(xiSol,1))
  }


  if(!is.null(options$multistart)){
    optRess<-NULL
    optRes<-list(value=-Inf,aux=NA)
    for(i in seq(options$multistart)){
      if(i==1){
        startingPoint<-colMeans(vertices_poly)
      }else{
        startingPoint<-lp(direction = "max",objective.in = runif(d-p),const.mat = Az,const.dir = "<",const.rhs = bz)$solution
      }
      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      #      optRess<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")
      optRess<-optim(par = startingPoint,fn = ff,gr = ffprime,method="BFGS")
      optRess$value= -optRess$value
      #      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
      if(optRess$value>=optRes$value){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    #    optRes<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")#,method = mm)
    optRes<-optim(par = startingPoint,fn = ff,gr = ffprime,method="BFGS")
    optRes$value = -optRes$value
  }

  optRes$solution = xiSol+Null(t(Phi))%*%optRes$par


  return(list(val=optRes$value,aux=optRes))
}


# getProfileInf function
#' @author Dario Azzimonti
#' @name getProfileInf_optim
#' @title Generic profile inf function computation with optim
#' @description Compute profile inf function for an arbitrary matrix \code{Phi} with with the L-BFGS-B algorithm of \link[stats]{optim}. Here the linear equality constraint is eliminated by using the Null space of \code{Phi}.
#' @param eta \eqn{p} dimensional point where the function is to be evaluated
#' @param Phi projection matrix of dimension \code{pxd}
#' @param f function to be optimized (takes a vector y of dimension d and returns a real number)
#' @param fprime derivative of f (same format, returning a \eqn{d} dimensional vector)
#' @param d dimension of the input for f
#' @param options a list containing the options to be passed to optim:
#' \itemize{
#' \item{\code{par:}}{contains the starting point (a point in dimension d)}
#' \item{\code{lower:} }{the lower bounds for the optimization domain (see optim for details)}
#' \item{\code{upper:}}{the upper bounds for the optimization domain (see optim for details)}
#' }
#' @return a real value corresponding to \eqn{min_{x \in D_Psi} f(x)}
#' @seealso \link{getProfileSup_optim}, \link{getProfileInf}, \link{plotMaxMin}
#' @export
getProfileInf_optim = function(eta,Phi,f,fprime,d,options=NULL){

  p = dim(Phi)[1]

  # put to zero small values
  Phi[which(abs(Phi)<.Machine$double.neg.eps)]<-0

  if(!is.null(options$trace))
    trace=options$trace

  # Find any solution of the linear system
  xiSol<-solutionLinConstraint(eta=eta,Phi=Phi,d = d,p=p)

  nullSpacePhi<-Null(t(Phi))

  # find equations for boundary conditions on z
  lower<-rep(0,d)
  upper<-rep(1,d)

  Az <- rbind(-nullSpacePhi,nullSpacePhi)
  bz <- rbind(xiSol-lower,-xiSol+upper)
  #  vertices_poly<-matrix(enumerate.vertices(A=Az,b=bz),ncol=d-p)
  rcdd_Hrep<-makeH(a1 = Az,b1 = bz)
  vertices_poly<-matrix(scdd(rcdd_Hrep)$output[,-c(1,2)],ncol=d-p)

  # generate point inside
  #  library(lpSolve)



  ff<-function(z){
    x_cand<-xiSol+nullSpacePhi%*%z
    g_cand<- sqrt(.Machine$double.xmax)
    if(all(x_cand>=rep(0,d))&all(x_cand<=rep(1,d))){
      return(f(x_cand))
    }else{
      return(g_cand)
    }
  }

  ffprime<-function(z){
    x_cand<-xiSol+nullSpacePhi%*%z
    if(all(Az%*%z-bz<=0)){
      return(crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))
    }else{
      ww<-Az%*%z-bz
      dirs<-which(ww>0)
      ww<-matrix(ww[dirs]/sum(ww[dirs]))
      return(-Az[dirs,]%*%ww)
    }
    return(crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))
  }

  if(!is.null(options$par)){
    #    startingPoint<-options$par
    #   if(Phi%*%startingPoint< eta)
    warning("Constraints check not implemented!")
  }else{
    #    startingPoint<-runif(min = lower,max = upper,d-p)
    startingPoint<-colMeans(vertices_poly)#lp(direction = "max",objective.in = runif(d-p),const.mat = Az,const.dir = "<=",const.rhs = bz)$solution
    #      matrix(c(xiSol,1))
  }


  if(!is.null(options$multistart)){
    #  cat("not implemented yet!")
    optRess<-NULL
    optRes<-list(value=Inf,aux=NA)
    for(i in seq(options$multistart)){
      #      startingPoint<-runif(d-p,min = lower,max=upper)
      if(i==1){
        startingPoint<-colMeans(vertices_poly)
      }else{
        startingPoint<-lp(direction = "max",objective.in = runif(d-p),const.mat = Az,const.dir = "<",const.rhs = bz)$solution
      }

      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      #      optRess<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")
      optRess<-optim(par = startingPoint,fn = ff,gr = ffprime,method="BFGS")
      if(optRess$value<=optRes$value){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    #    optRes<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")#,method = mm)
    optRes<-optim(par = startingPoint,fn = ff,gr = ffprime,method="BFGS")
  }


  optRes$solution = xiSol+Null(t(Phi))%*%optRes$par


  return(list(val=optRes$value,aux=optRes))

}

