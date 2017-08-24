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
#' @importFrom nloptr nloptr
#' @importFrom lhs maximinLHS
#' @importFrom splines bs
#' @importFrom graphics polygon
#' @importFrom methods is
#' @importFrom grDevices heat.colors terrain.colors
#' @importFrom stats as.formula
#' @importFrom RColorBrewer brewer.pal
#' @references
#'
#' Azzimonti, D., Bect, J., Chevalier, C., and Ginsbourger, D. (2016). Quantifying uncer- tainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D., Ginsbourger, D., Rohmer, J. and Idier, D. (2017+). Profile extrema for visualizing and quantifying uncertainties on excursion regions. Application to coastal flooding. Arxiv.
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Chevalier, C., Picheny, V., Ginsbourger, D. (2014). An efficient and user-friendly implementation of batch-sequential inversion strategies based on kriging. Computational Statistics & Data Analysis, 71: 1021-1034.
#'
#' Johnson, S. G.  The NLopt nonlinear-optimization package, http://ab-initio.mit.edu/nlopt
#'
#' Koenker, R. (2017). quantreg: Quantile Regression. R package version 5.33.
#'
#' Nocedal, J. and Wright, S. J. (2006). Numerical Optimization, second edition. Springer- Verlag, New York.
#'
#' Neuwirth, E. (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2.
#'
#' Roustant, O., Ginsbourger, D., Deville, Y. (2012). DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization. Journal of Statistical Software, 51(1): 1-55.
#'
NULL
