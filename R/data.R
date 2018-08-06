#' Coastal flooding as function of offshore conditions.
#'
#' A dataset containing the results of a numerical simulation conducted with
#' the MARS model (Lazure and Dumas, 2008) for coastal flooding.
#' The numerical model was adapted to the Boucholeurs area (French Atlantic coast),
#' close to La Rochelle, and validated with data from the 2010 Xynthia storm event.
#' See Azzimonti et al. (2017+) and Rohmer et al. (2018) for more details.
#'
#' The data frame contains 5 input variables: \code{Tide}, \code{Surge}, \code{phi}, \code{t-}, \code{t+} detailing
#' the offshore forcing conditions for the model. All input variables are normalized
#' in \eqn{[0,1]}. The response is \code{Area}, the area flooded in m^2.
#'
#'
#' @format A data frame with 200 rows and 6 variables:
#' \describe{
#'   \item{Tide}{High tide level in meters;}
#'   \item{Surge}{Surge height in meters;}
#'   \item{phi}{Phase difference between tide and surge;}
#'   \item{t-}{Duration of the increasing part of surge;}
#'   \item{t+}{Duration of the decreasing part of surge;}
#'   \item{Area}{Flooded area in m^2.}
#' }
#' @references
#'
#' Azzimonti, D., Ginsbourger, D., Rohmer, J. and Idier, D. (2017+). \emph{Profile extrema for visualizing and quantifying uncertainties on excursion regions. Application to coastal flooding.} arXiv:1710.00688.
#'
#' Rohmer, J., Idier, D., Paris, F., Pedreros, R., and Louisor, J. (2018). \emph{Casting light on forcing and breaching scenarios that lead to marine inundation: Combining numerical simulations with a random-forest classification approach.} Environmental Modelling & Software, 104:64-80.
#'
#' Rohmer, J. and Idier, D. (2012). \emph{A meta-modelling strategy to identify the critical offshore conditions for coastal flooding.} Natural Hazards and Earth System Sciences, 12(9):2943-2955.
#'
#' Lazure, P. and Dumas, F. (2008). \emph{An external-internal mode coupling for a 3D hydrodynamical model for applications at regional scale (MARS)}. Advances in Water Resources, 31:233-250.
#'
"coastal_flooding"
