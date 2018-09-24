#' WeMix: Package to Estimate Weighted Mixed-Effects Models. 
#' 
#' The WeMix package estimates mixed-effects models (also called multilevel models, 
#' mixed models, or hierarchical linear models) with survey weights. The likelihood function 
#' of such models with complex survey weights is not analytically calculable, so WeMix uses
#' numerical integration (Gauss-Hermite and adaptive Gauss-Hermite quadrature) to 
#' estimate mixed-effects models with survey weights at all levels of the model. 
#' 
#' This method allows users to analyze data that may have unequal selection probability at both 
#' the individual and group levels. Note that lme4 is the preferred way to estimate such 
#' models when there are no survey weights or weights only at the lowest level, and our 
#' estimation starts with parameters estimated in lme4. WeMix is intended for use in cases 
#' where there are weights at all levels,and is only for use with 
#' fully nested data. 
#' 
#' To start using WeMix, see the vignettes covering
#' the mathematical background of mixed-effects model estimation and use the
#' \code{mix} function to estimate models. Use 
#' \code{browseVignettes(package="WeMix")} to see the vignettes.
#' 
#' @docType package
#' @name WeMix-package
NULL
