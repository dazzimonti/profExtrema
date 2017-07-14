gradKm_dnewdata<-function (object, newdata, type, se.compute = TRUE, cov.compute = FALSE, 
          light.return = FALSE, bias.correct = FALSE,...) 
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
  der.formula<-deriv(object@trend.formula,namevec = colnames(object@X),function.arg = TRUE)
  trendDeriv<-function(x){
    result<-matrix(NA,ncol=ncol(x),nrow=nrow(x))
    for(i in seq(nrow(x)))
      result[i,]<-attr(do.call(der.formula,as.list(x[i,])),"gradient")
    return(cbind(rep(0,nrow(x)),result))
  }
  F.newdata <- trendDeriv(newdata) #model.matrix(trendDeriv, data = data.frame(newdata))
  y.predict.trend<-matrix(NA,nrow=nrow(newdata),ncol=ncol(newdata))
  y.predict.complement<-matrix(NA,nrow=nrow(newdata),ncol=ncol(newdata))
  for(i in seq(nrow(newdata))){
    if(length(beta)>1){
      y.predict.trend[i,] <- F.newdata[i,-1] * beta[-1]
    }else{
      y.predict.trend[i,] <- F.newdata[i,-1]
    }
    c.newdata <- covMat1Mat2(object@covariance, X1 = X, X2 = t(newdata[i,]), 
                             nugget.flag = object@covariance@nugget.flag)
    c.deriv.newdata<-covVector.dx(object = object@covariance,x = newdata[i,],X = X,c = c.newdata)
    Tinv.c.newdata <- backsolve(t(T), c.deriv.newdata, upper.tri = FALSE)
    y.predict.complement[i,] <- t(Tinv.c.newdata) %*% z
  }
  
  y.predict <- y.predict.trend + y.predict.complement
  output.list <- list()
  output.list$trend <- y.predict.trend
  output.list$mean <- y.predict
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
#      s2.predict.mat <- backsolve(t(T.M), t(F.newdata -
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
#       s2.predict.mat <- backsolve(t(T.M), t(F.newdata -
#                                               t(Tinv.c.newdata) %*% M), upper.tri = FALSE)
#       cond.cov <- cond.cov + crossprod(s2.predict.mat)
#       if (bias.correct)
#         cond.cov <- cond.cov * object@n/(object@n - object@p)
#     }
#     output.list$cov <- cond.cov
#   }
  return(output.list)
}