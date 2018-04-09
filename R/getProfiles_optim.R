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
# requires library(MASS)
getProfileSup_optim = function(eta,Phi,f,fprime,d,options=NULL){

  p = dim(Phi)[1]

  # put to zero small values
  Phi[which(abs(Phi)<.Machine$double.neg.eps)]<-0

  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace

  # Find any solution of the linear system
  #arbRest<-rep(0,d-p)#runif(d-p)
  #xiSol<-matrix(c(solve(Phi[1:p,1:p],eta-Phi[1:p,(p+1):d]%*%arbRest),arbRest))

  nullSpacePhi<-Null(t(Phi))

  # Find appropriate bounds for z
  # lower<- rep(0,d-p)
  # arbRest<-rep(1,d-p)#runif(d-p)
  # solvable<-Phi[1:p,1:p]
  # rest<-Phi[1:p,(p+1):d]
  # tempSol<-matrix(c(solve(solvable,eta-rest%*%arbRest),arbRest))
  # upper<-((tempSol-xiSol)/nullSpacePhi)[1,]

 lower <-rep(NA, d-p)
 upper <- rep(NA,d-p)
 allCombs<-combn(1:d,p)
 possRest<-numeric(0)
 xiSol=NULL
 indxx=NULL
 for(dd in seq(choose(d,p))){
   arbRest<-rep(0,d-p)#runif(d-p)
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

     if(is.null(indxx)){
       allNullSind<-combn(d,d-p)
       for(nn in seq(choose(d,d-p))){
         indxx<-allNullSind[,nn]
         if(kappa(nullSpacePhi[indxx,1:(d-p)])<1e16)
           break
       }
     }
     tempPossRest<-solve(nullSpacePhi[indxx,1:(d-p)],reOrderedSol[indxx,]-xiSol[indxx]) #((reOrderedSol-xiSol)/nullSpacePhi)[1,]
     possRest=rbind(possRest,tempPossRest)
   }
   arbRest<-rep(1,d-p)
   tempSol<-matrix(c(solve(solvable,eta-rest%*%arbRest),arbRest))
   reOrderedSol<- tempSol
   reOrderedSol[allCombs[,dd]]<-tempSol[1:p]
   reOrderedSol[setdiff(1:d,allCombs[,dd])]<-tempSol[(p+1):d]
   reOrderedSol<-round(reOrderedSol,digits = 12)
   if(all(reOrderedSol>=rep(0,d))&all(reOrderedSol<=rep(1,d))){
     if(is.null(xiSol))
       xiSol<-reOrderedSol

     if(is.null(indxx)){
       allNullSind<-combn(d,d-p)
       for(nn in seq(choose(d,d-p))){
         indxx<-allNullSind[,nn]
         if(kappa(nullSpacePhi[indxx,1:(d-p)])<1e16)
           break
       }
     }

     tempPossRest<-solve(nullSpacePhi[indxx,1:(d-p)],reOrderedSol[indxx,]-xiSol[indxx])#((reOrderedSol-xiSol)/nullSpacePhi)[1,]
     possRest=rbind(possRest,tempPossRest)
   }


 }
 lower=apply(possRest,2,min)
 upper=apply(possRest,2,max)



  # nloptr does minimizations so we need to change the sign
  ff<-function(z){return(-f(xiSol+nullSpacePhi%*%z))}
  ffprime<-function(z){return(-crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))}

  if(!is.null(options$par)){
#    startingPoint<-options$par
 #   if(Phi%*%startingPoint< eta)
      warning("Constraints check not implemented!")
  }else{
    startingPoint<-runif(min = lower,max = upper,d-p)
  }

#  if(d-p>1){
#    mm<-"BFGS"
#  }else{
#    mm<-"Brent"
#  }
  if(!is.null(options$multistart)){
    #  cat("not implemented yet!")
    optRess<-NULL
    optRes<-list(value=-Inf,aux=NA)
    for(i in seq(options$multistart)){
      startingPoint<-runif(d-p,min = lower,max=upper)
      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      optRess<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")
      optRess$value= -optRess$value
      #      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
      if(optRess$value>=optRes$value){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    optRes<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")#,method = mm)
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


  if(!is.null(options$lower))
    lower=options$lower
  if(!is.null(options$upper))
    upper=options$upper
  if(!is.null(options$trace))
    trace=options$trace

  # Find any solution of the linear system and appropriate bounds for z
  nullSpacePhi<-Null(t(Phi))


  lower <-rep(NA, d-p)
  upper <- rep(NA,d-p)
  allCombs<-combn(1:d,p)
  possRest<-numeric(0)
  xiSol=NULL
  indxx=NULL
  for(dd in seq(choose(d,p))){
    arbRest<-rep(0,d-p)#runif(d-p)
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

      if(is.null(indxx)){
        allNullSind<-combn(d,d-p)
        for(nn in seq(choose(d,d-p))){
          indxx<-allNullSind[,nn]
          if(kappa(nullSpacePhi[indxx,1:(d-p)])<1e16)
            break
        }
      }

      tempPossRest<-solve(nullSpacePhi[indxx,1:(d-p)],reOrderedSol[indxx,]-xiSol[indxx]) #((reOrderedSol-xiSol)/nullSpacePhi)[1,]
      possRest=rbind(possRest,tempPossRest)
    }
    arbRest<-rep(1,d-p)
    tempSol<-matrix(c(solve(solvable,eta-rest%*%arbRest),arbRest))
    reOrderedSol<- tempSol
    reOrderedSol[allCombs[,dd]]<-tempSol[1:p]
    reOrderedSol[setdiff(1:d,allCombs[,dd])]<-tempSol[(p+1):d]
    reOrderedSol<-round(reOrderedSol,digits = 12)
    if(all(reOrderedSol>=rep(0,d))&all(reOrderedSol<=rep(1,d))){
      if(is.null(xiSol))
        xiSol<-reOrderedSol

      if(is.null(indxx)){
        allNullSind<-combn(d,d-p)
        for(nn in seq(choose(d,d-p))){
          indxx<-allNullSind[,nn]
          if(kappa(nullSpacePhi[indxx,1:(d-p)])<1e16)
            break
        }
      }

      tempPossRest<-solve(nullSpacePhi[indxx,1:(d-p)],reOrderedSol[indxx,]-xiSol[indxx])#((reOrderedSol-xiSol)/nullSpacePhi)[1,]
      possRest=rbind(possRest,tempPossRest)
    }


  }
  lower=apply(possRest,2,min)
  upper=apply(possRest,2,max)



  # nloptr does minimizations so we need to change the sign
  ff<-function(z){return(f(xiSol+nullSpacePhi%*%z))}
  ffprime<-function(z){return(crossprod(fprime(xiSol+nullSpacePhi%*%z),nullSpacePhi))}

  if(!is.null(options$par)){
    #    startingPoint<-options$par
    #   if(Phi%*%startingPoint< eta)
    warning("Constraints check not implemented!")
  }else{
    startingPoint<-runif(min = lower,max = upper,d-p)
  }

  #  if(d-p>1){
  #    mm<-"BFGS"
  #  }else{
  #    mm<-"Brent"
  #  }
  if(!is.null(options$multistart)){
    #  cat("not implemented yet!")
    optRess<-NULL
    optRes<-list(value=Inf,aux=NA)
    for(i in seq(options$multistart)){
      startingPoint<-runif(d-p,min = lower,max=upper)
      if(!is.null(options$par) && i==1)
        startingPoint<-options$par
      optRess<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")
      #      optRess<-optim(par = startingPoint,fn = fxFixed,gr = fprimexFixed,control = list(fnscale=-1,trace=trace,pgtol=1e-4),method = methd,lower = lower,upper = upper)
      if(optRess$value<=optRes$value){
        optRes=optRess
      }else{
        break
      }
    }
  }else{
    optRes<-optim(par = startingPoint,fn = ff,gr = ffprime,lower=lower,upper=upper,method="L-BFGS-B")#,method = mm)
  }


  optRes$solution = xiSol+nullSpacePhi%*%optRes$par

  return(list(val=optRes$value,aux=optRes))

}

