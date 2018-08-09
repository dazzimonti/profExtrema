
<!-- README.md is generated from README.Rmd. Please edit that file -->
profExtrema
===========

Computation and plots of profile extrema functions. The package main functions are:

### Computation:

-   `coordinateProfiles`: Given a km objects computes the coordinate profile extrema function for the posterior mean and its quantiles.

-   `coordProf_UQ`: UQ part of coordinateProfiles.

-   `obliqueProfiles`: Given a km objects computes the profile extrema functions for a generic list of matrices Psi for the posterior mean and its quantiles.

-   `obliqueProf_UQ`: The UQ part of obliqueProfiles.

-   `getAllMaxMin`: computes coordinate profile extrema with full optimization for a deterministic function.

-   `approxMaxMin`: approximates coordinate profile extrema for a deterministic function.

-   `getProfileExtrema`: computes profile extrema given a list of matrices Psi for a deterministic function.

-   `approxProfileExtrema`: approximates profile extrema given a list of matrices Psi for a deterministic function.

### Plotting:

-   `plot_univariate_profiles_UQ`: plots for the results of coordProf\_UQ or obliqueProf\_UQ. Note that this function only works for univariate profiles.

-   `plotBivariateProfiles`: plots the bivariate maps results of a call to obliqueProfiles with a two dimensional projection matrix Psi.

-   `plotMaxMin`: simple plotting function for univariate profile extrema.

-   `plotOneBivProfile`: simple plotting function for bivariate profile extrema.

Note
----

This work was supported in part the Hasler Foundation, grant number 16065 and by the Swiss National Science Foundation, grant number 167199. The author warmly thanks David Ginsbourger, Jérémy Rohmer and Déborah Idier for fruitful discussions and accurate, thought provoking suggestions.

References
----------

Azzimonti, D., Bect, J., Chevalier, C., and Ginsbourger, D. (2016). Quantifying uncertainties on excursion sets under a Gaussian random field prior. SIAM/ASA Journal on Uncertainty Quantification, 4(1):850–874.

Azzimonti, D., Ginsbourger, D., Rohmer, J. and Idier, D. (2017+). Profile extrema for visualizing and quantifying uncertainties on excursion regions. Application to coastal flooding. arXiv:1710.00688.

Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.

Chevalier, C., Picheny, V., Ginsbourger, D. (2014). An efficient and user-friendly implementation of batch-sequential inversion strategies based on kriging. Computational Statistics & Data Analysis, 71: 1021-1034.

Johnson, S. G. The NLopt nonlinear-optimization package, <http://ab-initio.mit.edu/nlopt>

Koenker, R. (2017). quantreg: Quantile Regression. R package version 5.33.

Nocedal, J. and Wright, S. J. (2006). Numerical Optimization, second edition. Springer- Verlag, New York.

Neuwirth, E. (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2.

Roustant, O., Ginsbourger, D., Deville, Y. (2012). DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization. Journal of Statistical Software, 51(1): 1-55.
