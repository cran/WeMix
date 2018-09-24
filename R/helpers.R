# This helper funciton returns just the coefficient from the WeMix results object
#' @export
coef.WeMixResults <- function(object, ...) {
  object$coef
}

# This helper function prints the coefficient from the WeMix results object
#' @export
print.WeMixResults <- function(x, ...) {
  print(x$coef)
}

# This helper funciton creates a coefficient matrix from WeMix results
#' @export
summary.WeMixResults <- function(object, ...) {
  VC <- makeSandwich(object) # calculates the sandwhich estimator of variance 
  x0 <- object
  object$vc <- VC$vc #sets variances returned by model to variances calculated by sandwich estimator
  se <- VC$se[1:length(x0$coef)]
  object$coef <- cbind(Estimate=x0$coef, "Std. Error"=se, "t value"=x0$coef/se)
  rownames(object$coef) <- names(x0$coef)
  object$vars <- cbind(variance=x0$vars, "Std. Error"=VC$se[-(1:length(x0$coef))],"Std.Dev."=sqrt(x0$vars)) 
  class(object) <- "summaryWeMixResults"
  object
}

# This helper funciton creates a coefficient matrix from WeMix results
#' @export
#' @importFrom stats printCoefmat
print.summaryWeMixResults <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nVariance terms:\n")
  print(x$vars)
  cat("\nFixed Effects:\n")
  printCoefmat(x$coef)
  cat("\nlnl=",format(x$lnl,nsmall=2),"\n")
  cat("Intraclass Correlation=",format(x$ICC, nsmall=3, digits=3),"\n")
}


# This function calculates the sandwich estimator of the standard errors 
# standard error of estimators is calcuated using a sandwich estimator following Rabe-Hesketh and Skrondal 2006
#' @importFrom numDeriv jacobian
makeSandwich <- function(fittedModel) {

  par <- c(fittedModel$coef, fittedModel$vars)
  FisherInfInverse <- solve(fittedModel$hessian)
  
  L <- jacobian(fittedModel$lnlf, par, top=FALSE)
  nGroups <- nrow(L)
  nParams <- length(par)
  J <- matrix(0,ncol=nParams, nrow=nParams)
  for (i in 1:nGroups){
    J <- J + L[i,] %*% t(L[i,])
  }
  J <- (nGroups/(nGroups-1))*J
  SE <- FisherInfInverse%*%J%*%FisherInfInverse
  colnames(SE) <- rownames(SE)<-names(par)
  se <- sqrt(diag(SE)) 

  #return NA for obs with variance below threshold 
  se[which(log(fittedModel$vars) < -3.6) + length(fittedModel$coef) ] <- NA
  diag(SE) <- se^2 

  names(se) <- colnames(SE)
  return(list(VC=SE,se=se))
}
