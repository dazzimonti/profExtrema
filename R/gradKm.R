#' @author Dario Azzimonti
#' @name gradKm_dnewdata
#' @title Gradient of posterior mean and variance
#' @description Computes the gradient of the posterior mean and variance of the kriging model in \code{object} at the points \code{newdata}.
#' @param object a \link[DiceKriging]{km} object
#' @param newdata a vector, matrix or data frame containing the points where to perform predictions.
#' @param type a character corresponding to the type of kriging family (\code{"SK"} or \code{"UK"}).
#' @param se.compute an optional boolean indicating whether to compute the posterior variance or not. Default is TRUE.
#' @param bias.correct an optional boolean to correct bias in the UK variance. Default is FALSE.
#' @param light.return an optional boolean indicating whether to return additional variables. Default is FALSE.
#' @return Returns a list containing \itemize{
#' \item{\code{mean:}}{ the gradient of the posterior mean at \code{newdata}.}
#' \item{\code{trend:}}{ the gradient of the trend at \code{newdata}.}
#' \item{\code{s2:}}{ the gradient of the posterior variance at \code{newdata}.}
#' }
#' @export
gradKm_dnewdata<-function (object, newdata, type, se.compute = TRUE, light.return = FALSE, bias.correct = FALSE)
{
  nugget.flag <- object@covariance@nugget.flag
  X <- object@X
  y <- object@y
  newdata <- as.matrix(newdata)
  d.newdata <- ncol(newdata)
  if (!identical(d.newdata, object@d)) {
    stop("newdata must have the same numbers of columns than the experimental design")
  }
  if (!identical(colnames(newdata), colnames(X))) {
    colnames(newdata) <- colnames(X)
  }
  T <- object@T
  z <- object@z
  M <- object@M
  beta <- object@trend.coef

  F.newdata<-model.matrix(object@trend.formula, data = data.frame(newdata))

  # Pre process formula to remove things that are not in derivative tables
  newFormula<-as.formula(gsub(")","",gsub("I(","",object@trend.formula,fixed=TRUE),fixed=TRUE))

  der.formula<-deriv(newFormula,namevec = colnames(object@X),function.arg = TRUE)
  trendDeriv<-function(x){
    result<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
    for(i in seq(nrow(x)))
      result[i,]<-attr(do.call(der.formula,as.list(x[i,])),"gradient")
    return(cbind(rep(0,nrow(x)),result))
  }
  F.newdata.dx <- trendDeriv(newdata) #model.matrix(trendDeriv, data = data.frame(newdata))
  y.predict.trend<-matrix(NA,nrow=nrow(newdata),ncol=ncol(newdata))
  y.predict.complement<-matrix(NA,nrow=nrow(newdata),ncol=ncol(newdata))

  if(length(beta)>1){
    y.predict.trend <-sweep(F.newdata.dx,MARGIN = 2,beta,"*")[,-1]
  }else{
    y.predict.trend <- F.newdata.dx[,-1]
  }

  Tinv.c.newdata<-array(NA,dim = c(nrow(X),ncol(X),nrow(newdata)))
  c.deriv.newdata<-array(NA,dim = c(nrow(X),ncol(X),nrow(newdata)))
  c.newdata <- covMat1Mat2(object@covariance, X1 = X, X2 = newdata,
                           nugget.flag = object@covariance@nugget.flag)

  for(i in seq(nrow(newdata))){
    c.deriv.newdata[,,i]<-covVector.dx(object = object@covariance,x = newdata[i,],X = X,c = c.newdata[,i])
    Tinv.c.newdata[,,i] <- backsolve(t(T), c.deriv.newdata[,,i], upper.tri = FALSE)
    y.predict.complement[i,] <- crossprod(Tinv.c.newdata[,,i], z)
  }

  y.predict <- y.predict.trend + y.predict.complement
  output.list <- list()
  output.list$trend <- y.predict.trend
  output.list$mean <- y.predict

  if(se.compute){
    if (!is(object@covariance, "covUser")) {
      total.sd2 <- 0
    }
    else {
      stop("User covariance kernel derivatives not implemented")
    }

    s2.predict.1<-matrix(NA,nrow=nrow(newdata),ncol=ncol(newdata))
    for(i in seq(nrow(newdata))){
      s2.predict.1[i,]<-crossprod(Tinv.c.newdata[,,i],backsolve(t(T),c.newdata[,i],upper.tri = FALSE))
    }
    if (type == "SK") {
      s2.predict <- total.sd2 - 2*s2.predict.1
    }else if (type == "UK") {
      T.M <- chol(t(M) %*% M)
      s2.predict.mat <- matrix(NA,nrow=nrow(newdata),ncol=ncol(newdata))
      for(i in seq(nrow(newdata))){
        s2.predict.mat[i,]<-crossprod(backsolve(t(T.M), t(F.newdata[i,] -  backsolve(t(T),c.newdata[,i])%*% M), upper.tri = FALSE),
                                      backsolve(t(T.M), t(F.newdata.dx[i,-1] - crossprod(Tinv.c.newdata[,,i], M)), upper.tri = FALSE))
      }
      s2.predict <- total.sd2 - 2*s2.predict.1 + 2*s2.predict.mat
      if (bias.correct)
        s2.predict <- s2.predict * object@n/(object@n - object@p)
    }
#    lower95 <- y.predict - q95 * sqrt(s2.predict)
#    upper95 <- y.predict + q95 * sqrt(s2.predict)
    output.list$s2 <- s2.predict
#    output.list$lower95 <- lower95
#    output.list$upper95 <- upper95
  }
#  if (!light.return) {
#    output.list$c <- c.newdata
#    output.list$Tinv.c <- Tinv.c.newdata
# }
#  if ((se.compute) || (cov.compute)) {
#    if (!is(object@covariance, "covUser")) {
#      total.sd2 <- object@covariance@sd2
#    }
#    else {
#      m <- nrow(newdata)
#      total.sd2 <- rep(NA, m)
#      for (i in 1:m) {
#        total.sd2[i] <- object@covariance@kernel(newdata[i,
#                                                        ], newdata[i, ])
#     }
#   }
#    if (object@covariance@nugget.flag) {
#      total.sd2 <- total.sd2 + object@covariance@nugget
#    }
#  }
#  if (se.compute) {
#    s2.predict.1 <- apply(Tinv.c.newdata, 2, crossprod)
#    if (type == "SK") {
#      s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
#      s2.predict <- as.numeric(s2.predict)
#      q95 <- qnorm(0.975)
#    }
#    else if (type == "UK") {
#      T.M <- chol(t(M) %*% M)
#      s2.predict.mat <- backsolve(t(T.M), t(F.newdata.dx -
#                                              t(Tinv.c.newdata) %*% M), upper.tri = FALSE)
#      s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
#      s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2,
#                         0)
#      s2.predict <- as.numeric(s2.predict)
#      if (bias.correct)
#        s2.predict <- s2.predict * object@n/(object@n -
#                                                object@p)
#       q95 <- qt(0.975, object@n - object@p)
#     }
#     lower95 <- y.predict - q95 * sqrt(s2.predict)
#     upper95 <- y.predict + q95 * sqrt(s2.predict)
#     output.list$sd <- sqrt(s2.predict)
#     output.list$lower95 <- lower95
#     output.list$upper95 <- upper95
#   }
#   if (cov.compute) {
#     C.newdata <- covMatrix(object@covariance, newdata)[[1]]
#     cond.cov <- C.newdata - crossprod(Tinv.c.newdata)
#     if (type == "UK") {
#       T.M <- chol(t(M) %*% M)
#       s2.predict.mat <- backsolve(t(T.M), t(F.newdata.dx -
#                                               t(Tinv.c.newdata) %*% M), upper.tri = FALSE)
#       cond.cov <- cond.cov + crossprod(s2.predict.mat)
#       if (bias.correct)
#         cond.cov <- cond.cov * object@n/(object@n - object@p)
#     }
#     output.list$cov <- cond.cov
#   }
  return(output.list)
}
