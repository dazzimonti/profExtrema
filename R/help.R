#' @title profExtrema package
#' @description Computation and plots of profile extrema functions. The package main functions are: \describe{
#'    \item{\strong{Computation:}}{ \itemize{
#'     \item \code{\link{coordinateProfiles}}: Given a \link[DiceKriging]{km} objects computes the coordinate profile extrema function for the posterior mean and its quantiles.
#'     \item \code{\link{coordProf_UQ}}: UQ part of \code{coordinateProfiles}.
#'     \item \code{\link{obliqueProfiles}}: Given a \link[DiceKriging]{km} objects computes the profile extrema functions for a generic list of matrices Psi for the posterior mean and its quantiles.
#'     \item \code{\link{obliqueProf_UQ}}: The UQ part of \code{obliqueProfiles}.
#'     \item \code{\link{getAllMaxMin}}: computes coordinate profile extrema with full optimization for a deterministic function.
#'     \item \code{\link{approxMaxMin}}: approximates coordinate profile extrema for a deterministic function.
#'     \item \code{\link{getProfileExtrema}}: computes profile extrema given a list of matrices Psi for a deterministic function.
#'     \item \code{\link{approxProfileExtrema}}: approximates profile extrema given a list of matrices Psi for a deterministic function.
#'    } }
#'    \item{\strong{Plotting:}}{ \itemize{
#'     \item \code{\link{plot_univariate_profiles_UQ}}: plots for the results of \code{coordProf_UQ} or \code{obliqueProf_UQ}. Note that this function only works for univariate profiles.
#'     \item \code{\link{plotBivariateProfiles}}: plots the bivariate maps results of a call to \code{obliqueProfiles} with a two dimensional projection matrix Psi.
#'     \item \code{\link{plotMaxMin}}: simple plotting function for univariate profile extrema.
#'     \item \code{\link{plotOneBivProfile}}: simple plotting function for bivariate profile extrema.
#'    }  }
#'
#' }
#' @details Package: profExtrema \cr
#' Type: Package \cr
#' Version: 0.2.1 \cr
#' Date: 2020-03-20
#'
#' @author Dario Azzimonti (dario.azzimonti@@gmail.com) .
#' @docType package
#' @name profExtrema
#' @import microbenchmark
#' @importFrom DiceKriging predict.km covMat1Mat2 covVector.dx covMatrix simulate km
#' @importFrom KrigInv integration_design
#' @importFrom pGPx krig_weight_GPsimu optim_dist_measure grad_kweights
#' @importFrom quantreg rq predict.rq
#' @importFrom grDevices cairo_pdf dev.off gray.colors adjustcolor pdf
#' @importFrom graphics abline contour image legend lines par plot points hist axis
#' @importFrom stats deriv model.matrix optim predict rnorm quantile runif smooth.spline
#' @importFrom utils combn
#' @importFrom utils modifyList
#' @importFrom MASS Null
#' @importFrom rcdd makeH scdd lpcdd d2q q2d
#' @importFrom lhs maximinLHS
#' @importFrom splines bs
#' @importFrom graphics polygon layout
#' @importFrom methods is
#' @importFrom grDevices heat.colors terrain.colors
#' @importFrom stats as.formula
#' @importFrom RColorBrewer brewer.pal
#' @note This work was supported in part the Hasler Foundation, grant number 16065 and by the Swiss National Science Foundation, grant number 167199. The author warmly thanks David Ginsbourger, Jérémy Rohmer and Déborah Idier for fruitful discussions and accurate, thought provoking suggestions.
#' @references
# @importFrom lpSolve lp
#'
#' Azzimonti, D., Bect, J., Chevalier, C., and Ginsbourger, D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.
#'
#' Azzimonti, D., Ginsbourger, D., Rohmer, J. and Idier, D. (2017+). Profile extrema for visualizing and quantifying uncertainties on excursion regions. Application to coastal flooding. arXiv:1710.00688.
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
