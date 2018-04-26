# getProfileSup function with smarter algorithm!
#' @author Dario Azzimonti
#' @name solutionLinConstraint
#' @title Obtain generic solution to linear constraints
#' @description Computes an arbitrary solution \eqn{x^*} that satisfies \code{Phi x = eta}
#' @param eta \eqn{p} dimensional point where the function is to be evaluated
#' @param Phi \eqn{p x d} matrix describing the projection matrix
#' @param low lower limits of the box constraints on the input (vector of dimension d)
#' @param upp upper limits of the box constraints on the input (vector of dimension d)
#' @param pars a vector of dimension d-p that can be selected by the user to guide the solution. When NULL (suggested) it is chosen by the function.
#' @return a vector of dimensions \code{d x 1} containing the solution \eqn{x^*}.
# requires library(MASS) and lpSolve
solutionLinConstraint<- function(eta,Phi,low=NULL,upp=NULL,pars=NULL){
  d<-ncol(Phi)
  p<-nrow(Phi)

  if(is.null(low)){
    low=rep(0,d)
  }

  if(is.null(upp)){
    upp=rep(1,d)
  }

  kk<-0
  xiSol<-NULL
  allCombs<-combn(1:d,p)

  while(is.null(xiSol)){
    kk<-kk+1
    for(dd in seq(choose(d,p))){
      if(is.null(pars)){
        if(kk==1){
          arbRest<-rep(0,d-p)
        }else if(kk==2){
          arbRest<-rep(1,d-p)
        }else{
          arbRest<-runif(d-p)
        }
      }else{
        arbRest=pars
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

      if(all(reOrderedSol>=low)&all(reOrderedSol<=upp)){
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

  # find equations for boundary conditions on z
  ## TODO: add option to set box-constraints
  lower<-rep(0,d)
  upper<-rep(1,d)

  # Verify if the feasible set is not empty
  Hlarge<-makeH(a1 = rbind(diag(-1,nrow = d,ncol=d),diag(1,nrow = d,ncol=d)),
                b1 = rbind(matrix(lower,nrow = d),matrix(upper,nrow = d)),
                a2 = Phi,b2=matrix(eta))

  Hlarge<-d2q(Hlarge)
  solLP<-lpcdd(hrep = Hlarge,objgrd = d2q(rep(0,d)))

  if(solLP$solution.type=="Optimal"){
    # Find any solution of the linear system
    xiSol<-matrix(q2d(solLP$primal.solution),nrow=d) #solutionLinConstraint(eta=eta,Phi=Phi,low = lower,upp=upper)
  }else{
    warning("Cannot proceed with these constraints. A simple LP optimization problem does not have an optimal solution")
    optRes<-list(value=NA,solLP=solLP,message="The feasible set is either empty (solLP$solution.type==Inconsistent) or the problem is dual inconsistent. See solLP for the results of a LP problem with those contraints and objective function 0.")
    return(list(val=optRes$value,aux=optRes))
  }
  # rm(Hlarge,solLP)
  # gc()


  nullSpacePhi<-Null(t(Phi))



  Az <- rbind(-nullSpacePhi,nullSpacePhi)
  bz <- rbind(xiSol-lower,-xiSol+upper)
  #  vertices_poly<-matrix(enumerate.vertices(A=Az,b=bz),ncol=d-p)
  rcdd_Hrep0<-makeH(a1 = Az,b1 = bz)
  rcdd_Hrep<-d2q(rcdd_Hrep0)
  # rm(rcdd_Hrep0)
  # gc()
  vertices_poly<-matrix(q2d(scdd(rcdd_Hrep)$output[,-c(1,2)]),ncol=d-p)


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
      return(crossprod(ww,-Az[dirs,]))
      #return(-Az[dirs,]%*%ww)
    }
    return(-crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))
  }

  if(!is.null(options$par)){
    #    startingPoint<-options$par
    #   if(Phi%*%startingPoint< eta)
    warning("Constraints check not implemented!")
  }else{
    #    startingPoint<-runif(min = lower,max = upper,d-p)
    startingPoint<-colMeans(vertices_poly)
    #      matrix(c(xiSol,1))
  }

  if(dim(vertices_poly)[1]==1){
    # The feasible region for the subproblem is one point, the starting point,
    # so it doesn't make sense to try different starting points
    options$multistart=NULL
  }

  if(!is.null(options$multistart)){
    optRess<-NULL
    optRes<-list(value=-Inf,aux=NA)
    for(i in seq(options$multistart)){
      if(i==1){
        startingPoint<-colMeans(vertices_poly)
      }else if(i==2){
        m_lp<-lpcdd(rcdd_Hrep,objgrd = d2q(rep(0,d-p)))
        if(m_lp$solution.type=="Optimal")
          startingPoint<-q2d(m_lp$primal.solution)
      }else if(i==3){
        m_lp<-lpcdd(rcdd_Hrep,objgrd = d2q(rep(1,d-p)),minimize = FALSE)
        if(m_lp$solution.type=="Optimal")
          startingPoint<-q2d(m_lp$primal.solution)
      }else{
        m_lp<-lpcdd(rcdd_Hrep,objgrd = d2q(runif(d-p)))
        if(m_lp$solution.type=="Optimal")
          startingPoint<-q2d(m_lp$primal.solution)
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

  rm(rcdd_Hrep0,Hlarge,solLP)
#  gc()

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


  # find equations for boundary conditions on z
  ## TODO: add option to set box-constraints
  lower<-rep(0,d)
  upper<-rep(1,d)

  # Verify if the feasible set is not empty
  Hlarge<-makeH(a1 = rbind(diag(-1,nrow = d,ncol=d),diag(1,nrow = d,ncol=d)),
                b1 = rbind(matrix(lower,nrow = d),matrix(upper,nrow = d)),
                a2 = Phi,b2=matrix(eta))
  Hlarge<-d2q(Hlarge)
  solLP<-lpcdd(hrep = Hlarge,objgrd = d2q(rep(0,d)))

  if(solLP$solution.type=="Optimal"){
    # Find any solution of the linear system
    xiSol<-matrix(q2d(solLP$primal.solution),nrow=d) #solutionLinConstraint(eta=eta,Phi=Phi,low = lower,upp=upper)
  }else{
    warning("Cannot proceed with these constraints. A simple LP optimization problem does not have an optimal solution")
    optRes<-list(value=NA,solLP=solLP,message="The feasible set is either empty (solLP$solution.type==Inconsistent) or the problem is dual inconsistent. See solLP for the results of a LP problem with those contraints and objective function 0.")
    return(list(val=optRes$value,aux=optRes))
  }
  # rm(Hlarge,solLP)
  # gc()


  nullSpacePhi<-Null(t(Phi))


  Az <- rbind(-nullSpacePhi,nullSpacePhi)
  bz <- rbind(xiSol-lower,-xiSol+upper)
  #  vertices_poly<-matrix(enumerate.vertices(A=Az,b=bz),ncol=d-p)
  rcdd_Hrep0<-makeH(a1 = Az,b1 = bz)
  rcdd_Hrep<-d2q(rcdd_Hrep0)
  vertices_poly<-matrix(q2d(scdd(rcdd_Hrep)$output[,-c(1,2)]),ncol=d-p)


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
      return(crossprod(ww,-Az[dirs,]))
      #return(-Az[dirs,]%*%ww)
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

  if(dim(vertices_poly)[1]==1){
    # The feasible region for the subproblem is one point, the starting point,
    # so it doesn't make sense to try different starting points
    options$multistart=NULL
  }


  if(!is.null(options$multistart)){
    #  cat("not implemented yet!")
    optRess<-NULL
    optRes<-list(value=Inf,aux=NA)
    for(i in seq(options$multistart)){
      if(i==1){
        startingPoint<-colMeans(vertices_poly)
      }else if(i==2){
        m_lp<-lpcdd(rcdd_Hrep,objgrd = d2q(rep(0,d-p)))
        if(m_lp$solution.type=="Optimal")
          startingPoint<-q2d(m_lp$primal.solution)
      }else if(i==3){
        m_lp<-lpcdd(rcdd_Hrep,objgrd = d2q(rep(1,d-p)),minimize = FALSE)
        if(m_lp$solution.type=="Optimal")
          startingPoint<-q2d(m_lp$primal.solution)
      }else{
        m_lp<-lpcdd(rcdd_Hrep,objgrd = d2q(runif(d-p)))
        if(m_lp$solution.type=="Optimal")
          startingPoint<-q2d(m_lp$primal.solution)
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

  rm(rcdd_Hrep0,Hlarge,solLP)
#  gc()

  return(list(val=optRes$value,aux=optRes))

}

