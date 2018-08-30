#' Coastal flooding as function of offshore forcing conditions.
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
#'   \item{Surge}{Surge peak amplitude in meters;}
#'   \item{phi}{Phase difference between high tide and surge peak;}
#'   \item{t-}{Duration of the increasing part of the surge temporal signal (assumed to be triangular);}
#'   \item{t+}{Duration of the decreasing part of the surge temporal signal (assumed to be triangular);}
#'   \item{Area}{Flooded area in m^2.}
#' }
#'
#' @examples
#' # Define inputs
#' inputs<-data.frame(coastal_flooding[,-6])
#' colnames(inputs)<-colnames(coastal_flooding[,-6])
#' colnames(inputs)[4:5]<-c("tPlus","tMinus")
#'
#' # put response in areaFlooded variable
#' areaFlooded<-data.frame(coastal_flooding[,6])
#' colnames(areaFlooded)<-colnames(coastal_flooding)[6]
#' response = sqrt(areaFlooded)
#'
#' \donttest{
#' model <- km(formula=~Tide+Surge+I(phi^2)+tMinus+tPlus,
#'             design = inputs,response = response,covtype="matern3_2")
#' # Fix threshold
#' threshold<-sqrt(c(1.2e6,1.9e6,3.1e6,6.5e6))
#'
#' # use the coordinateProfile function
#' ## set up plot options
#' options_plots <- list(save=FALSE, folderPlots = "./" ,
#'                       titleProf = "Coordinate profiles",
#'                       title2d = "Posterior mean",qq_fill=TRUE)
#' # set up full profiles options
#' options_full<-list(multistart=15,heavyReturn=TRUE)
#' # set up approximation options
#' d <- model@@d
#' init_des<-lhs::maximinLHS(5*d , d )
#' options_approx<- list(multistart=2,heavyReturn=TRUE, initDesign=init_des,
#'                       fullDesignSize=100, smoother="quantSpline")
#' # set up simulation options
#' ## reduce nsims and batchsize for faster/less accurate UQ
#' nsims=200
#' opts_sims<-list(algorithm="B", lower=rep(0,d ),
#'                 upper=rep(1,d ), batchsize=150,
#'                 optimcontrol=list(method="genoud", pop.size=100,print.level=0),
#'                 integcontrol = list(distrib="sobol",n.points=2000),nsim=nsims)
#'
#' opts_sims$integration.param <- integration_design(opts_sims$integcontrol,
#'                                                   d , opts_sims$lower,
#'                                                   opts_sims$upper,
#'                                                   model,threshold)
#' opts_sims$integration.param$alpha <- 0.5
#' # run the coordinate profile extrema on the mean
#' CF_CoordProf_mean<- coordinateProfiles(object = model, threshold = threshold,
#'                                        uq_computations = FALSE, options_approx = options_approx,
#'                                        plot_level=3, plot_options= options_plots, return_level=3,
#'                                        options_full=options_full)
#'
#' # run UQ computations
#' CF_CoordProf_UQ<- coordinateProfiles(object = CF_CoordProf_mean, threshold = threshold,
#'                                      uq_computations = TRUE, options_approx = options_approx,
#'                                      plot_level=3, plot_options= options_plots, return_level=3,
#'                                      options_sims=opts_sims,options_full=options_full,
#'                                      options_bound = list(beta=0.024,alpha=0.05))
#' }
#'
#' @references
#'
#' Azzimonti, D., Ginsbourger, D., Rohmer, J. and Idier, D. (2017+). \emph{Profile extrema for visualizing and quantifying uncertainties on excursion regions. Application to coastal flooding.} arXiv:1710.00688.
#'
#' Rohmer, J., Idier, D., Paris, F., Pedreros, R., and Louisor, J. (2018). \emph{Casting light on forcing and breaching scenarios that lead to marine inundation: Combining numerical simulations with a random-forest classification approach.} Environmental Modelling & Software, 104:64-80.
#'
#' Lazure, P. and Dumas, F. (2008). \emph{An external-internal mode coupling for a 3D hydrodynamical model for applications at regional scale (MARS)}. Advances in Water Resources, 31:233-250.
#'
"coastal_flooding"
