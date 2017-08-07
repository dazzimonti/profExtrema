#' @title profExtrema package
#' @description Computation of profile extrema functions. The package main functions are: \itemize{
#'     \item \code{\link{coordinateProfiles}}: Given a \link[DiceKriging]{km} objects computes the profile extrema function for the posterior mean an its quantiles.
#'     \item \code{\link{coordProf_UQ}}: UQ part of \code{coordinateProfiles}.
#'     \item \code{\link{getAllMaxMin}}: central function to compute the coordinate profile extrema with full optimization.
#'     \item \code{\link{approxMaxMin}}: central function to compute the coordinate profile extrema with approximations.
#' }
#' @details Package: profExtrema \cr
#' Type: Package \cr
#' Version: 0.1.0 \cr
#' Date: 2017-07-14
#'
#' @author Dario Azzimonti (dario.azzimonti@@gmail.com) .
#' @docType package
#' @name profExtrema
#' @import microbenchmark
#' @importFrom DiceKriging predict.km covMat1Mat2 covVector.dx covMatrix simulate
#' @importFrom KrigInv integration_design
#' @importFrom GPsims krig_weight_GPsimu optim_dist_measure grad_kweights simulate_km
#' @importFrom quantreg rq predict.rq
#' @importFrom grDevices cairo_pdf dev.off gray.colors adjustcolor
#' @importFrom graphics abline contour image legend lines par plot points
#' @importFrom stats deriv model.matrix optim predict rnorm quantile runif smooth.spline
# @importFrom nloptr nloptr
#' @importFrom lhs maximinLHS
#' @importFrom splines bs
#' @importFrom graphics polygon
#' @importFrom methods is
#' @importFrom grDevices heat.colors terrain.colors
#' @importFrom stats as.formula
#' @importFrom RColorBrewer brewer.pal
#' @references TODO
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
NULL
