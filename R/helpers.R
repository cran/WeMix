# This helper funciton returns just the coefficient from the WeMix results object.
#' @method coef WeMixResults
#' @export
coef.WeMixResults <- function(object, ...) {
  object$coef
}

# This helper function prints the coefficient from the WeMix results object.
#' @export
print.WeMixResults <- function(x, ...) {
  print(x$coef)
}

# This helper funciton creates a coefficient matrix from WeMix results.
# and packages variance results in a clearer way. 
#' @importFrom stats cov2cor model.frame
#' @export
summary.WeMixResults <- function(object, ...) {
  x0 <- object
  
  #for adaptive quad models 
  if(object$is_adaptive){
    VC <- makeSandwich(object) # calculates the sandwhich estimator of variance 
    object$vc <- VC$VC #sets variances returned by model to variances calculated by sandwich estimator
    se <- VC$se[1:length(x0$coef)]
    re_se <- VC$se[-(1:length(x0$coef))] #not used currently 
    x0$varVC <- list(NULL,VC$VC[names(x0$vars),names(x0$vars)]) #This is hard coded for 2 levels, with a NULL because the are no refs at level 1
    # build the var mat output, first adding SEvcov
    varDF <- x0$varDF
    varDF$SEvcov <- NA
    names(re_se) <- gsub(":",".",names(re_se),fixed=TRUE) # change : to .
    for(i in 1:nrow(varDF)) {
      if(varDF$fullGroup[i] %in% names(re_se)) {
        varDF$SEvcov[i] <- re_se[varDF$fullGroup[i]] 
      }
    }
  } else{
    se <- object$SE
    re_se <- rep(0,length(x0$vars))
    # build the var mat output from this
    varDF <- x0$varDF
  }
  object$coef <- as.matrix(data.frame(Estimate=x0$coef, "Std. Error"=se, "t value"=x0$coef/se))
  # make.names breaks/fixes colnames with spaces, un-break/-fix them.
  colnames(object$coef)[2:3] <- c("Std. Error", "t value")
  rownames(object$coef) <- names(x0$coef)
  # build the var mat output
  varsmat <- varDF[is.na(x0$varDF$var2),c("level","grp", "var1","vcov","SEvcov")]
  varsmat$st <- sqrt(varsmat$vcov)
  # add variance estimates
  colnames(varsmat) <- c("Level", "Group", "Name", "Variance", "Std. Error", "Std.Dev.")
  for(li in 2:x0$levels) {
    vc <- as.matrix(x0$varVC[[li]])
    cr <- cov2cor(vc)
    if(ncol(vc)>1) {
      for(i in 2:ncol(cr)) {
        for(j in 1:(i-1)){
          varsmat[varsmat$Level==li & varsmat$Name==rownames(cr)[i],paste0("Corr",j)] <- cr[i,j]
        }
      }
    }
  }

  object$varsmat <- varsmat
  #this is important for compatibility  with mixed.sdf
  object$vars  <- varsmat[,4:6]
  rownames(object$vars)  <- names(x0$vars)
  
  class(object) <- "summaryWeMixResults"
  return(object)
}

# This helper function prints WeMix Wald Test Resuls.
#' @export
print.WeMixWaldTest <- function(x, ...) {
  cat("Wald Test Statistic\n")
  cat(x$W)
  cat("\n\np-value\n")
  cat(round(x$p,4))
  cat("\n\nDegrees of Freedom\n")
  cat(x$df)
  cat("\n\nNull Hypothesis\n")
  if (class(x$H0)=="matrix"){
    print(apply(x$H0,FUN=round,digits=4,MARGIN=c(1:length(dim(x$H0)))))
  } else {
    print(round(x$H0,4))
  }
  cat("\n\nAlternative Hypothesis\n")
  if (class(x$HA)=="matrix"){
    print(apply(x$HA,FUN=round,digits=4,MARGIN=c(1:length(dim(x$H0)))))
  } else {
    print(round(x$HA,4))
  }
}

# This helper function creates a coefficient matrix from WeMix results.
#' @export
#' @importFrom stats printCoefmat
print.summaryWeMixResults <- function(x, digits = max(3, getOption("digits") - 3), nsmall=2, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nVariance terms:\n")
  varsmat <- x$varsmat
  varsmat <- varsmat[order(varsmat$Level, decreasing=TRUE),]
  i <- 1
  while(paste0("Corr",i) %in% colnames(varsmat)) {
    # round correlations to two digits
    varsmat[[paste0("Corr",i)]] <- as.character(round(varsmat[[paste0("Corr",i)]],2))
    i <- i + 1
  }
  print(varsmat, na.print="", row.names=FALSE, digits=digits, nsmall=nsmall)
  cat("Groups:\n")
  groupSum <- varsmat[!duplicated(varsmat$Level), c("Level", "Group")]
  if(any(groupSum$Level==1)) {
    groupSum$Group[groupSum$Level==1] <- "Obs"
  } else {
    newrow <- groupSum[1,]
    newrow$Level <- 1
    newrow$Group <- "Obs"
    groupSum <- rbind(groupSum, newrow)
  }
  groupSum$"n size" <- rev(x$ngroups)
  for(i in 1:length(x$wgtStats)) {
    groupSum$"mean wgt"[groupSum$Level == i] <- x$wgtStats[[i]]$mean
    groupSum$"sum wgt"[groupSum$Level == i] <- x$wgtStats[[i]]$sum
  }
  print(groupSum, na.print="", row.names=FALSE, digits=digits, nsmall=nsmall)
  cat("\nFixed Effects:\n")
  printCoefmat(x$coef, digits=digits, ...)
  cat("\nlnl=",format(x$lnl, nsmall=nsmall),"\n")
  if(length(x$ICC) > 0) {
    cat("Intraclass Correlation=",format(x$ICC, nsmall=nsmall, digits=digits),"\n")
  }
}

#' Mixed Model Wald Tests
#' @description
#' This function calculates the Wald test for either fixed effects or variance parameters.
#' @param fittedModel a  model of class \code{WeMixResults} that is the result of a call to \code{\link{mix}} 
#' @param type a string, one of "beta" (to test the fixed effects) or "Lambda" (to test the variance-covariance parameters for the random effects)
#' @param coefs a vector containing the names of the coefficients to test. For \code{type="beta"} these must be
#'              the variable names exactly as they appear in the fixed effects table of the summary. For \code{type="Lambda"}
#'              these must be the names exactly as they appear in the theta element of the fitted model. 
#' @param hypothesis the hypothesized values of beta or Lambda. If \code{NA} (the default) 0 will be used. 
#' @details
#' By default this function tests against the null hypothesis that all coefficients are zero.
#' To identify which coefficients to test use the name exactly as it appears in 
#' the summary of the object.
#' @return Object of class \code{WeMixWaldTest}. 
#' This is a list with the following elements: 
#' \item{wald}{the value of the test statistic.}
#' \item{p}{the p-value for the test statistic. Based on the probabilty of the test statistic
#' under the chi-squared distribution.}
#' \item{df}{degrees of freedom used to calculate p-value.}
#' \item{H0}{The vector (for a test of beta) or matrix (for tests of Lambda) containing the 
#' null hypothesis for the test.}
#' \item{HA}{The vector (for a test of beta) or matrix (for tests of Lambda) containing the 
#' alternative hypothesis for the test (i.e. the values calculated by the fitted model being 
#' tested.)} 
#' @examples 
#' \dontrun{
#' library(lme4) #to use the example data 
#' sleepstudyU <- sleepstudy
#' sleepstudyU$weight1L1 <- 1
#' sleepstudyU$weight1L2 <- 1
#' wm0 <- mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, 
#'             weights=c("weight1L1", "weight1L2"))
#' wm1 <- mix(Reaction ~ Days + (1 +Days|Subject), data=sleepstudyU, 
#'            weights=c("weight1L1", "weight1L2"))
#' 
#' waldTest(wm0, type="beta")  #test all betas 
#' #test only beta for days 
#' waldTest(wm0, type="beta", coefs="Days")  
#' #test only beta for intercept against hypothesis that it is 1
#' waldTest(wm0, type="beta", coefs="(Intercept)", hypothesis=c(1))  
#' 
#' waldTest(wm1,type="Lambda")  #test all values of Lambda
#'  #test only some Lambdas.The names are the same as names(wm1$theta)
#' waldTest(wm1,type="Lambda", coefs="Subject.(Intercept)") 
#' #specify  test values
#' waldTest(wm1,type="Lambda", coefs="Subject.(Intercept)", hypothesis=c(1))  
#' }
#' @importFrom stats pchisq 
#' @export
waldTest <- function(fittedModel, type=c("beta", "Lambda") , coefs=NA, hypothesis=NA) {
  type <- match.arg(type,c("beta", "Lambda"))
  if (type == "Lambda"){
    # names of the coefs 
    if (all(is.na(coefs))) {
      coef_names <- names(fittedModel$theta)
    } else {
      coef_names <- coefs
    } 
    
    #number of parameters 
    q <- length(coef_names)
    
    #total coefficients 
    p <- length(fittedModel$theta)
    
    # this block sets up R, the hypothesis matrix 
    # it has a row for each parameter to be tested (q) 
    # it has a column for each parameter that exists in the model
    R <-  matrix(0, nrow=q, ncol=p)
    rownames(R) <- coef_names
    colnames(R) <- names(fittedModel$theta)
    
    # Fill in one for each  each coefficient to  be tested in the relevant row
    for (coef in coef_names){
      if (any(!coef_names %in% names(fittedModel$theta))){
        stop("Names of coefficients to test must be the same as names of theta in the fitted model.")
      }
      R[coef, coef] <- 1
    }
    
    #this block sets up the r vector of hypothesized values for theta
    if (all(is.na(hypothesis))){
      r <- rep(0,q)
    } else {
      if (length(hypothesis) != q){
        stop("Length of hypothesized values must be same as number of coefficients to test.")
      }
      r <- hypothesis
    }
    names(r) <- coef_names
    
    # get variance covariance matrix 
    V_hat <- fittedModel$var_theta

    #estimates 
    ests <- fittedModel$theta
    
    #create var-cov matrix for null and alternative hypothesis 
    variances <- fittedModel$varDF[!is.na(fittedModel$varDF$var1) & is.na(fittedModel$varDF$var2),"fullGroup"]
    h0  <- matrix(0,nrow=length(variances),ncol=length(variances))
    rownames(h0) <- colnames(h0)  <- variances
    ha <- h0
   
    #fill with theta values (for null hypothesis) - but only those related to the parameters  we are test
    #handle positioning covariance
    for (i  in 1:length(r)){
      #decompose name into parts
      var <- r[i]
      parts <- unlist(strsplit(names(var),".",fixed=T))
      if (length(parts)>2){  #These are covariances 
        group <- parts[1]
        v1 <-paste0(group,".",parts[2])
        v2 <-paste0(group,".",parts[3])
      }  else { #these are variances 
        group <- parts[1]
        v1 <-paste0(group,".",parts[2])
        v2 <-paste0(group,".",parts[2])
      }
      #put value alue into var-covar matrix
      h0[v1,v2] <- var
      h0[v2,v1] <- var
      
      ha[v1,v2] <- fittedModel$theta[names(var)]
      ha[v2,v1] <- fittedModel$theta[names(var)]
    }
  } else { # end if (type == "Lambda")
    # wald test for beta values
    # names of the coefs 
    if  (all(is.na(coefs))) {
      coef_names <- names(fittedModel$coef)
    } else {
      coef_names <- coefs
    } 
    
    #number of parameters 
    q <- length(coef_names)
    
    #total coefficients 
    p <- length(fittedModel$coef)
    
    # this block sets up  R the hypothesis matrix 
    # it has a row for each parameter to be tested (q) 
    # it has a column for each parameter that exists in the model
    R <-  matrix(0,nrow=q,ncol=p)
    rownames(R) <- coef_names
    colnames(R) <- names(fittedModel$coef)
    
    # Fill in one for each  each coefficient to  be tested in the relevant row
    for (coef in coef_names){
      if (any(!coef_names %in% names(fittedModel$coef))){stop("Names of coefficients to test must be the same as names of beta in the fitted model.")}
      R[coef,coef] <- 1
    }
    
    #this block sets up the r vector of hypothesized values for theta
    if (all(is.na(hypothesis))){
      r <- rep(0,q)
    } else {
      if (length(hypothesis) != q){stop("Length of hypothesized values must be same as number of coefficients to test.")}
      r <- hypothesis
    }
    
    # Covariance matrix of beta
    if(fittedModel$is_adaptive){
      # this is a glm
      VC <- makeSandwich(fittedModel) # calculates the sandwhich estimator of variance 
      V_hat <- VC$VC[1:(q+1), 1:(q+1)]
    } else {
      V_hat <- fittedModel$cov_mat
    }
   
    #estimates 
    ests <- fittedModel$coef
    
    #Make hypothesis vectors to return
    h0 <- ests[coef_names]
    ha <- r 
    names(ha) <- coef_names
  } #end else for if (type == "beta")
  
  #Wald test using the following equation
  #W = (Rb - r)^T (RVR^T)^-1 (Rb - r)
  W = as.numeric(t(R%*%ests - r)  %*% solve(R  %*% V_hat  %*% t(R)) %*% (R%*%ests - r))
  ch = 1 - pchisq(W, df=q)
  res <- list("Wald"=W,"p"=ch ,"df"=q,"H0"= h0,"HA"=ha)
  class(res) <- "WeMixWaldTest"
  return(res)
}  

# This function calculates robust standard errors using the sandwich estimator of the standard errors  
# following Rabe-Hesketh & Skrondal 2006. For linear models, robust standard errors are returned from the main estimation function. 
# This function is only used for post-hoc estimation for non-linear models. 
#' @importFrom numDeriv jacobian
makeSandwich <- function(fittedModel) {
  #make same estimator based on 

  par <- c(fittedModel$coef, fittedModel$vars)
  FisherInfInverse <- solve(fittedModel$invHessian)
  
  
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

# TRUE if x has no zeros in it. False if it does have a zero.
# useful for finding columns with non-zero entries in lapply
nozero <- function(x) {
  any(x != 0)
}
