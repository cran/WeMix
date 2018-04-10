#' Survey Weighted Mixed-Effects Models 
#' @param formula  a formula object in the style of lme4 that creates the model.
#' @param data  a data frame containing the raw data for the model.
#' @param weights a character vector of weight variables found in data frame.
#' @param nQuad an integer  number of quadrature point to evaluate on. See notes for guidelines. 
#' @param run  logical, true runs the model, false provides partial output for debug or testing.
#' @param keepAdapting logical, set to TRUE when the adaptive quadrature should adapt after every Newton step. Defaults to TRUE. 
#' False should be used rarely for faster (but less accurate) results.
#' @param verbose boolean, default False, set to TRUE to print results of intermediate steps 
#' @param acc0 integer, the precision of mpfr, default 120 
#' @param start numeric vector representing the point at which the model should start optimization; takes the shape of c(coef, vars) 
#' from results (see help). 
#' @param fast logical, use  c++ function for faster result. Defaults to True. 
#' @description
#' Implements a survey weighted mixed-effects model using the provided formula. 
#' @details
#' Uses adaptive quadrature following the method in Statas's GLAMMM. For additional details, see main vignette 
#' "WeightedMixedModels.pdf"  which provides extensive examples as well as a description of the mathematical 
#' basis of the estimation procedure. The main specification also shows comparisons to model specifications in other
#' common software. 
#' 
#' Notes: 
#' \itemize{
#' \item Standard errors of random effect variances are estimated by the Sandwich Method, see main vingette for details. 
#' \item To see the function that is maximized in the estimation of this model see the section on "Model fitting" in the main vingette.
#' \item When all weights above the individual level are one, this is similar to a lmer and you should use lme4 because it is much faster.
#' \item Starting coefficients are not provided they are estimated using lme4. 
#' \item When the variance of a random effect is very low (<.1) we don't estimate it because very low variances create problems with  numerical evaluation.
#'  In these cases, consider re estimating without that RE.
#'  \item The model is estimated by maximum likelihood estimation, restricted maximum liklihood (REML) is not available 
#' \item To chose number of quadrature points a balance is needed between accuracy and speed- estimation time increases quadratically 
#' with the number of points chosen. Additionaly an odd number of points is traditionaly used. We recommend starting at 13 and increasing 
#' or decreasing as needed. 
#' }
#' @importFrom lme4 getME lmer 
#' @importFrom stats dnorm aggregate terms
#' @importFrom numDeriv genD hessian grad
#' @useDynLib WeMix, .registration = TRUE
#' @importFrom Rcpp evalCpp 
#' @return object of class WeMixResults
#' This is a list with objects: 
#' lnlf - function, the likelihood function 
#' lnl - numeric, the logliklihood of the model 
#' coef - numeric vector, the estimated coefficients of the model 
#' vars- numeric vector, the variances
#' call - the original call used 
#' levels - integer, the number of levels in the model 
#' ICC - numeric, the Intraclass Correlation Coefficient 
#' CMEAN  - function the conditional mode function that can be called with par and omega to get the conditional mode of the likelihood function
#' CMODE - function the conditional mean function that can be called with par and omega to get the conditional mean of the likelihood function
#' Hessian - the second derivative of the likelihood function 
#' @examples 
#' library(WeMix)
#' library(EdSurvey)
#' naep <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package="NAEPprimer"))
#' naepData <- getData(naep,varnames=c("lep","origwt","smsrswt","scrpsu","dsex","mrpcm1"))[1:200,]
#' me1 <- mix(mrpcm1 ~  lep+ dsex + (1|scrpsu), data=naepData, 
#'    weights=c("origwt", "smsrswt"), nQuad=7, verbose=FALSE, fast=TRUE,run=TRUE)
#' @author Paul Bailey, Claire Kelley, and Trang Nguyen 
#' @export
mix <- function(formula, data, weights, nQuad=13L, run=TRUE, verbose=TRUE, acc0=120, keepAdapting=FALSE, start=NULL, fast=FALSE) {
  #############################################
  #                   Outline:                #   
  #     1) data preparation and reshaping     #
  #    2) Identify integration parameters     #
  #     3) Maximum Likelihood Estimation      #
  #            4) Post Estimation             #
  #############################################
  
  #############################################
  #              1) Data Preparation          #
  #############################################
  call <- match.call()
  # check class and some requirements of arguments
  if(!inherits(formula, "formula")) stop(paste0("the argument ", sQuote("formula"), " must be a formula."))
  if(!inherits(data, "data.frame")) stop(paste0("the argument ", sQuote("data"), " must be a data.frame."))  
  if(nQuad <= 0) stop(paste0(sQuote("nQuad"), " must be a positive integer."))
  if(!inherits(run, "logical")) stop(paste0("the argument ", sQuote("run"), " must be a logical."))
  if(!inherits(verbose, "logical")) stop(paste0("the argument ", sQuote("verbose"), " must be a logical."))
  if(!inherits(weights, "character")) stop(paste0("The argument ", sQuote("weights"), " must be a character vector of weight column names in ", sQuote("data"), "."))
  if(any(!weights %in% colnames(data))) stop(paste0("The argument ", sQuote("weights"), " must specify valid columns in ", sQuote("data"), "."))
  if(acc0 <= 0) stop(paste0(sQuote("acc0"), " must be a positive integer."))
  
  #set up initial values 
  adapter <- "MAP" #initial function evaluation through MAP, BLUE estimator can be used post hoc
  weights0 <- weights # keep track of incomming weights column names
  # correctly format acc0
  acc0 <- round(acc0)
  nQuad <- round(nQuad) #nquad must be an integer
  
  # Manually calculate group names (ie level 2 variables) from the formula
  # 
  group_name <- c()
  termLabels <- grep("\\|", attr(terms(formula),"term.labels"), value = TRUE)
  group_name <- unique(trimws(gsub(".*\\|","", termLabels)))
  group_name <- unique(unlist(strsplit(group_name,"/")))
  
  # reorder data by groups (in order to make Z matrix obtained from lme easier to work with)
  for (i in length(group_name):1) {
    if (group_name[i] %in% names(data)) {
      data <- data[order(data[,group_name[i]]),]
    }
  }
  # remove row names so that resorted order is used in lme model 
  row.names(data) <- NULL
  
  # run lmer to get a ballpark fit and model structure information
  if(verbose) {
    cat("Using lmer to get an approximate (unweighted) estimate and model structure.\n")
  }
  lme <- lmer(formula, data, REML=FALSE)
  #Get the Z (random effects) matrix from LME 
  model_matrix <- getME(lme,"mmList")
  #Get names of group (level 2) varaibles for random effects
  #group_name <- unique(names(getME(lme, "Tp")[-1]))
  
  
  z_groups <- names(model_matrix) #get full names random effects whcih include both variable and group level
  
  # store the full sample weights in wgts0
  wgts0 <- data[,weights]
  
  # create columns for any interactions coming from the formula
  # this will be mainly used in 3 level models and is included for forward compatability 
  need_vars <- group_name[!group_name %in% names(data)] #find which group names are not in data set (ie composite groups)
  names(data) <- gsub(":","_",names(data)) #any original names have : stripped out so as not to confuse with composite groups
  if (length(need_vars)>0){
    for (i in 1:length(need_vars)){
      col_name <- need_vars[i]
      vars<-gsub(":","_",strsplit(need_vars[i],":")[[1]]) #also remove : here so that we match new col names
      data[,col_name]<-apply(data[,vars],1,function(x) {paste(x,collapse=":")}) #paste all level of interaction variables together
    }
  }
  
  # prepare Z matrices for each level 
  Z <- list(NULL)
  # Z is de-duplicated. ZFull is always the size of the outcome vector
  ZFull <- list(NULL)
  n_rows <- nrow(data)
  for (i in 1:length(group_name)){
    z_to_merge <- grepl(paste0("[|]\\W",group_name[i],"$"), z_groups)
    Z_i <- matrix( unlist(model_matrix[z_to_merge], use.names=FALSE), nrow=n_rows)
    ZFull <- c(ZFull, list(Z_i))
    if(i > 1) {
      Z_i <- Z_i[!duplicated(data[,group_name[i-1]]),,drop=FALSE]
    }
    Z <- c(Z, list(Z_i))
  }
  
  # find number of random effects classes
  nz <- list(0) # there are not REs at level 1
  for(i in 1:length(Z)) {
    if(!is.null(Z[[i]])) {
      nz[[i]] <- ncol(Z[[i]])
    }
  }
  
  #calculate n levels from Z and warn if not the same dimension as weights 
  levels <- length(Z)
  if(length(weights) != levels) {
    stop(paste0("the argument ", sQuote("weights"), " must be a list of column names with length equal to levels."))  
  }
  
  
  # transform weights into a list of  dataframes where each data frame has
  # columns representing unit weight at that level and index of group names 
  weights <- list()
  for(i in 1:length(nz)) {
    df <- data.frame(w=wgts0[,i],stringsAsFactors=FALSE)
    # add the index for this level, when there is one
    if(i < length(nz)) {
      df$indexp1 <- data[,group_name[i]]
    }
    if(i > 1) {
      # for levels >1 data frame indexed by group name and values represent weights 
      df$index <- data[,group_name[i-1]]
      df <- df[!duplicated(df$index),]
    }
    weights[[i]] <- df
  }

  #get the y variable name from the formula 
  y_label <- as.character(formula[[2]])
  
  # get lmer fixed effects and covariance terms
  k <- length(lmeb <- getME(lme, "fixef"))
  parlme <- c(lmeb)
  lmesummary <- summary(lme)
  
  # find number of unique groups 
  ngrp <- lmesummary$ngrps
  if(length(unique(ngrp)) != length(ngrp)) {
    # if non-nested groups are present (ie not all obs at lower levels are part of upper level groups) 
    # then there will be non unique entries in ngrp 
    stop("The model does not appear to have nested groups.")
  }

  # set up variance and coefficients 
  lmeVarDF <- data.frame(lmesummary$varcor)
  parlme <- c(parlme, lmeVarDF$vcov)
  
  # use initial values for coefficients if they are provied, otherwise default to lme4
  if(is.null(start)) {
    est0 <- parlme
  } else {
    if(length(start) != length(parlme)) {
      stop(paste0("Expecting argument ", sQuote("start"), " to have ", length(est0), " elements, found ", length (start), " elements."))
    }
    est0 <- start
    names(est0) <- names(parlme)
  } 
  
  # prepare variance data frame
  # rename groups in var-cov matrix to be all distinct 
  ind <- 1
  while(sum(grepl(paste0("\\.", ind, "$"), lmeVarDF$grp)) > 0) {
    lmeVarDF$grp <- sub(paste0("\\.", ind, "$"), "", lmeVarDF$grp)
    ind <- ind + 1
  }
  lmeVarDF$sdcor <- NULL # remove extraneous column
  lmeVarDF$ngrp <- NA
  lmeVarDF$grp <- gsub(".",":",lmeVarDF$grp,fixed=TRUE)
  
  # add column showing how many elements are in each group
  for(vari in 1:nrow(lmeVarDF)) {
    if(lmeVarDF$grp[vari] == "Residual") {
      lmeVarDF$ngrp[vari] <- nrow(data)
    } else {
      lmeVarDF$ngrp[vari] <- ngrp[names(ngrp) == lmeVarDF$grp[vari]]
    }
  }
  
  # add column indicating which model level each group belongs to
  ngrp2 <- rev(sort(unique(lmeVarDF$ngrp)))
  for(grpi in 1:length(ngrp2)) {
    lmeVarDF$level[ngrp2[grpi] == lmeVarDF$ngrp] <- grpi
  }
  
  # add names to the variance terms of the paramter vector
  names(est0)[-(1:k)] <- lmeVarDF$grp
  
  
  # use helper funtion to create a covariance matrix from the data frame with variance and covariance information
  covarianceConstructor <- covMat2Cov(lmeVarDF)
  # C is the realization of the covariance constructor
  C <- covarianceConstructor(est0[-(1:k)])
  
  # these are the realized y/X vector/matrix
  y <- data[,c(y_label)]
  X <- getME(lme,"X")
  
  ##########################################################
  #           2) Identify integration parameter            #
  ##########################################################
  
  if(verbose) {
    cat("Identifying initial integration locations (MAP) estimates for random effects.\n")
  }

  # setup these methods that, theoretically, can be used to find the 
  # adaptive integration points. For now, only MAP works.
  # Note: MAP finds the maximum of the likelihood surface.
  # it is not, properly, a MAP
  MAP0 <- MAP(groups=data[,group_name,drop=FALSE], y=y, X=X, levels=levels,
                  Z=Z, ZFull=ZFull, weights=weights, k=k,
                  qp=gauss.quad(nQuad, "hermite"),
                  covariance_constructor=covarianceConstructor, verbose=verbose,
                  nlmevar=nrow(lmeVarDF)-1, nz=nz, acc=acc0)
  # notes: BLUE finds the expected value of the likelihood surface
  # it is not, properly, a BLUE
  BLUE0 <- BLUE(groups=data[,group_name,drop=FALSE], y=y, X=X, levels=levels,
                Z=Z, ZFull=ZFull, weights=weights, k=k,
                qp=gauss.quad(nQuad, "hermite"),
                covariance_constructor=covarianceConstructor, verbose=verbose,
                nlmevar=nrow(lmeVarDF)-1, nz=nz, acc=acc0)
  
  # form omega0
  # omega0 is a list of matricies where each list element represents a level
  # each matrix column represents the value of b for all units at the level 
  bvec <- getME(lme, "b") #get b values from lme results
  bvecCuts <- getME(lme, "Gp") #find out indices corresponding to different levels
  blist <- vector("list", levels) # NULL for level 1
  
  #startLoc helps transform the z vector into a z matrix, starting from the beginning ( index of 1)
  startLoc <- 1
  
  # number of rows that the Z matrix has at level (l-1), skipped 1 here
  # this is used exclusively to form bmat
  # the issue is that bvecCuts doesn't have cuts for factor levels
  # the n_rows_z code fixes thatS
  comps <- names(getME(lme,"cnms")) #has same n as bvec and tells you group and var
  n_rows_z <- list()
  for (i in 1:length(comps)){
    n_rows_z[i] <- lmeVarDF[lmeVarDF$grp == comps[i],"ngrp"][1]
  }
  
  #forming matrixes of b values 
  blist <- vector("list", levels) # NULl for level 1
  for(cuti in 2:length(bvecCuts)) {
    bmat <- matrix(bvec[startLoc:bvecCuts[cuti]], nrow=n_rows_z[[cuti-1]])
    li <- unique(lmeVarDF$level[lmeVarDF$ngrp==nrow(bmat)])
    blist[[li]] <- cbind(blist[[li]], bmat)
    startLoc <- bvecCuts[cuti] + 1
  }
  
  omega0 <- blist
  
  # make MAP estimate with b values and parameter estimates 
  a0 <- MAP0(omega0=omega0, par0=est0)

  # add determinant to scale z values
  # created the adapted quadrature locaiton by moving original points based on curvature Q
  # follows equation of Hartzel et al, pg 87 (Reference in main vignette)
  zScale <- lapply(a0$Qi0, function(Qi0i) {
               if(is.null(Qi0i)) {
                 return(NULL)
               }
               df <- data.frame(detQ=sapply(Qi0i,det))
               # add group labels
               for(i in 1:length(group_name)) {
                 if(length(unique(data[,group_name[i]])) == nrow(df)) {
                   df[,group_name[i]] <- unique(data[,group_name[i]])
                   attr(df,"groups") <- c(attr(df, "groups"), group_name[i])
                 }  
               }
               df
             })
  
  index <- data.frame(data[,c(group_name)])
  names(index) <- group_name
  
  #add column with determinant of Q to the Z infomration (join by index which is group name)
  for(wi in 2:length(weights)) {
    Zgrps <- attr(zScale[[wi]], "groups")
    weights[[wi]] <- merge(weights[[wi]],zScale[[wi]][,c(Zgrps, "detQ")],by.x="index", by.y=Zgrps)
  }
  
  # first guess
  est <- est0

  # use lnl function to optimize
  qp <- gauss.quad(nQuad,"hermite")
  # the covariance constructor maps the variance terms that are less than 1
  # to exp(x-1) to avoid negative variance terms.
  # fn0 does not use this map and can be easily used by the user
  # fn0R does use that map and is intended for use by the optimizer
  # these calls are otherwise identical.

  fn0 <- param.lnl.quad(y=y,
                        X=X,
                        levels=levels,
                        Z=Z,
                        ZFull=ZFull,
                        Qi=a0$Qi,
                        QiFull=a0$QiFull,
                        omega=a0$omega,
                        omegaFull=a0$omegaFull,
                        W=weights,
                        k=k,
                        qp=qp,
                        cConstructor=covarianceConstructor,
                        acc0=acc0,
                        fast=fast,
                        mappedDefault=FALSE)
  
  fn0R <- param.lnl.quad(y=y,
                         X=X,
                         levels=levels,
                         Z=Z,
                         ZFull=ZFull,
                         Qi=a0$Qi,
                         QiFull=a0$QiFull,
                         omega=a0$omega,
                         omegaFull=a0$omegaFull,
                         W=weights,
                         k=k,
                         qp=qp,
                         cConstructor=covarianceConstructor,
                         acc0=acc0,
                         fast=fast,
                         mappedDefault=TRUE)

  if(!run) {
    # if run option is false the function terminates here, returning the components used in the optimization
    # without performing the optimization
    return(list(lnlf=fn0R, parlme=parlme, omega0=a0$omega0, lme=lme, adapt=a0, weights=weights))
  }
  
  #############################################
  #     3) Maximum Likelihood estimation      #
  #############################################
  
  d1 <- rep(Inf, length(est)) #d1 represent first derivatives of estimates 
  oldlnl <- fn0(est) # evaluate the likelihood with the estimated coefficients before adaptive quadrature 
  
  a00 <- a0 # start from previous estimates (MAP estimates) 
    
  if(verbose) {
    cat("Starting Newton steps\n")
  }
  
  #select just variances that are > -3 to aviod continued adaptation below 0 
  covs_and_vars <- est[-(1:k)]
  vars <- covs_and_vars[which(is.na(lmeVarDF$var2))]
  not_0_vars <- which(vars>-3)+k
  est[-(1:k)]  <- ifelse(est[-(1:k)]< -4.6,-4.6,est[-(1:k)])  #reset variance that get below threshold
  # v is the full Newton step vector. Here set to a stub value
  v <- d1

  # Newton steps continue until step size is under threshold, excluding variances < -3
  # the denominator prevents relative changes that are indistinguishable from 0 from being
  # included.
  while(max(abs(v[c(1:k,not_0_vars)]/pmax(abs(est[c(1:k,not_0_vars)]),1e-5))) > 1E-5) { 
    
    gradAndHessian <- getHessianAndGrad(fn0, est) # calls helper function to get both first and second derivative 
    d1 <- gradAndHessian$grad
    d2 <- gradAndHessian$hessian
    fact <- 1 # scaling factor by which this step is multiplied 
    v <- solve(d2) %*% d1 # the Newton step
    if(verbose) {
      cat("lnl:",oldlnl, " max (relative) derivative=", max(abs(v[c(1:k,not_0_vars)]/pmax(abs(est[c(1:k,not_0_vars)]),1e-5))), " ")
      cat("\ncurrent solution, gradient, and Newton step:\n")
      prnt <- cbind(oldEstimate=est, firstDeriv=d1, proposedNewtonEstimate=est -1*v)
      rownames(prnt) <- names(est0)
      colnames(prnt) <- c("previous Est", "firstDeriv", "Newton Step")
      print(prnt)
    }
    
    # create new estimate by stepping in the direction idicated by gradient and scaled by fact
    newest <- est - fact * v
    newlnl <- fn0(newest)
    stp <- 0
    
    # make sure Newton step is a step improves lnl
    while(newlnl < oldlnl) {
      stp <- stp + 1
      if(verbose) {
        cat("halving step\n")
      }
      fact <- fact/2
      if(stp > 5 & fact > 0) {
        if(verbose) {
          cat("reversing\n")
        }
        ##reverse direction if more than 5 steps have been taken, avoid newtons method getting stuck 
        fact <- -1
        stp  <- 0
      }
      if (stp>10) {
        warning("Newton's method gives a bad search direction, stopping optimization. Results may not have converged.")
        fact <- 0
        keepAdapting <- FALSE
        d1 <- d1 * 0
        v <- v * 0
        oldlnl <- -Inf
      }
      newest <- est - fact * v
      newlnl <- fn0(newest)      
    } # end while(newlnl < oldlnl)
    if(verbose) {
      cat("\n")
    }
    # update current estimate with the step we decided on
    est <- est - fact * v
    # update the lnl
    oldlnl <- newlnl

    # now, worry about threshold
    if(sum(est[-(1:k)]< -3.59) > 0) {
      est[-(1:k)]  <- ifelse(est[-(1:k)]< -3.59,-3.59,est[-(1:k)])  #reset variance that get below threshold
      oldlnl <- fn0(est)      
    }
    #adapts new quadrature points if parameter keepAdapting is true  
    if(keepAdapting) {
      if(verbose) {
        cat("adapting\n")
      }
      # adapter is never BLUE now
      if(adapter == "BLUE") {
        a0 <- BLUE0(omega0=a0$omega0, par0=est0, Qi0=a0$Qi0)
      } else {
        # update the b and Q parameters
        a0 <- MAP0(omega0=a0$omega0,
                   par0=est,
                   verb=FALSE)
      }
      
      #recacluate scaled Z values based on new Q (Curvature)
      zScale <- lapply(a0$Qi0, function(Qi0i) {
                   if(is.null(Qi0i)) {
                     return(NULL)
                   }
                   df <- data.frame(detQ=sapply(Qi0i,det))
                   # add group labels
                   for(i in 1:length(group_name)) {
                     if(length(unique(data[,group_name[i]])) == nrow(df)) {
                       df[,group_name[i]] <- unique(data[,group_name[i]])
                       attr(df,"groups") <- c(attr(df, "groups"), group_name[i])
                     }  
                   }
                   df
                 })
      
      # move new values of detQ on to dataframes in weights list 
      for(wi in 2:length(weights)) {
        weights[[wi]]$detQ <- NULL
        Zgrps <- attr(zScale[[wi]], "groups")
        weights[[wi]] <- merge(weights[[wi]],zScale[[wi]][,c(Zgrps, "detQ")],by.x="index", by.y=Zgrps)
      }

      # update function for unser and function for optimization
      fn0 <- param.lnl.quad(y=y,
                            X=X,
                            levels=levels,
                            Z=Z,
                            ZFull=ZFull,
                            Qi=a0$Qi,
                            QiFull=a0$QiFull,
                            omega=a0$omega,
                            omegaFull=a0$omegaFull,
                            W=weights,
                            k=k,
                            qp=qp,
                            cConstructor=covarianceConstructor,
                            acc0=acc0,
                            fast=fast,
                            mappedDefault=FALSE)
      
      fn0R <- param.lnl.quad(y=y,
                             X=X,
                             levels=levels,
                             Z=Z,
                             ZFull=ZFull,
                             Qi=a0$Qi,
                             QiFull=a0$QiFull,
                             omega=a0$omega,
                             omegaFull=a0$omegaFull,
                             W=weights,
                             k=k,
                             qp=qp,
                             cConstructor=covarianceConstructor,
                             acc0=acc0,
                             fast=fast,
                             mappedDefault=TRUE)
      
      # Process of adapting stops when there the relative difference between chosen points is below threshold 
      if(max(abs(a00$omega0[[2]] - a0$omega0[[2]])/pmax(abs(a0$omega0[[2]]),1E-10)) < 1E-2) {
        if(verbose) {
          cat("done adapting: the MAP is not changing sufficnelty.\n")
        }
        keepAdapting <- FALSE
      }
      if(keepAdapting & max(abs(d1)) <= 1E-3) {
        if(verbose) {
          cat("done adapting: close to solution.\n")
        }
        keepAdapting <- FALSE
      }
      a00 <- a0
      # the lnl function may have shifted enough that we need to update the current lnl
      # so update it based on these new quad points
      oldlnl <- fn0(est)
    } # end if(keepAdapting)
    #end of the Newton's method while loop
    
    #re updates the index of vars that are greater than -3 based on new est values 
    covs_and_vars <- est[-(1:k)]
    vars <- covs_and_vars[which(is.na(lmeVarDF$var2))] #select only vars not covars 
    not_0_vars <- which(vars > -3) + k
  } # end of while(max(abs(d1/pmax(abs(est),1e-5))) > 1E-5)

  #############################################
  #            4) Post Estimation             #
  #############################################
 
  hessian <- d2
  
  # For each "random effect" calculate the maximum of the likelihood surface (MAP)
  # and the expected value of the likelihood surface (BLUE)
  MAP <- MAP0(omega0=a0$omega0, par0=est, verb=FALSE)$omega0
  BLUE <- BLUE0(omega0=a0$omega0, par0=est,Qi0=a0$Qi0,adapt=FALSE,verb=FALSE)
  
  # make est numeric
  est <- as.numeric(est)
  
  # add final check to check whether fast method agrees with the slow (pure R, high precision)
  if (fast) {
    lnl <- fn0(est)
    if (abs(lnl - fn0(est,fast0 = FALSE)) > 0.001 * abs(lnl)) {
      warning("The likelihood function may be inaccurate, try fitting with fast=FALSE")
    }
  }
  # make sure the names agree with lmer
  names(est) <- names(parlme)
  
  # fix vars that are less than one to be mapped
  covs_and_vars <- est[-(1:k)]
  vars <- covs_and_vars[which(is.na(lmeVarDF$var2))] #select only vars not covars 
  need_fix_vars <- which(vars < 1)
  covs_and_vars[need_fix_vars] <- exp(covs_and_vars[need_fix_vars] - 1)
  
  #if any of variances  got re mapped, re calculate hessian with new position 
  if (length(need_fix_vars)>0){
    hessian <- d2(fn0,c(est[1:k],covs_and_vars))
  }

  #calculate interclass corelation (ICC) 
  vars <- covs_and_vars
  names(vars) <- gsub(":NA","",paste(lmeVarDF$grp,lmeVarDF$var1,lmeVarDF$var2,sep=":"))

  var_between <- sum(vars[which(!is.na(lmeVarDF$var1) & is.na(lmeVarDF$var2))])
  var_within <- vars[which(lmeVarDF$grp=="Residual")]
  ICC <- var_between/(var_between+var_within)
  
  # it is possible for the lnl to have changed slightly, so update it to avoid confusion
  
  res <- list(lnlf=fn0R, lnl=fn0(est), coef=est[1:k], vars=vars,
              call=call, levels=levels, CMEAN=MAP,ICC=ICC, CMODE=BLUE,
              hessian=hessian)
  class(res) <- "WeMixResults"
  return(res)
}


# finds expected value of the "random effects" likelihood surfaces
# this function cannot find these without input estimates; it cannot be used for adapting.
# @author Paul Bailey
BLUE <- function(groups, y, X, levels, Z, ZFull, weights0, k, qp, covariance_constructor, verbose, nlmevar, nz, acc, fast=FALSE) {

  # must define one of Qi or Qi0. Defining Qi is faster
  #this is the function that is returned when BLUE is called 
  function(omega0, par0, Qi=NULL, Qi0=NULL, verb=verbose, adapt=TRUE) {
    
    # form the Q matrixes from weights
    weights <- weights0 #start with original weights 
    if(is.null(Qi)) {
      Qi <- list(NULL)
      QiFull <- list(NULL)
      for( oi in 2:length(omega0)) {
        map <- groups[,oi-1] # decrement index by one because the first groups data frame doesnt have a column for level 1
        umap <- unique(map)
        nzi <- ncol(Z[[oi]]) # find number of z variables at this level 
        Qi[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(weights[[oi-1]]))
        for(i in 1:nrow(weights[[oi-1]])) {
           #shape Qi, a list with one element per model level
          #each element of Qi has one row for each random effect (Z) and one column for each observation in the level below-random effect combination
          #values are being moved over from Qi0 which is the MAP estimates at the group level
          Qi[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
        }
        QiFull[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
        for(i in 1:nrow(X)) {
          #Also  full version which as one row for each Z and one column for each obs-random effect combination
          QiFull[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
        }
      } # ends  for( oi in 2:length(omega0))
    } # ends if(is.null(Qi)

    # function that creates dataframe to be used to scale Z value by the scaling factor calculated as the determinate of the curvature Q
    zScale <- lapply(Qi0, function(Qi0i) {
                 if(is.null(Qi0i)) {
                   return(NULL)
                 }
                 df <- data.frame(detQ=sapply(Qi0i,det))
                 # add group labels 
                 for(i in 1:ncol(groups)) {
                   if(length(unique(groups[,i])) == nrow(df)) {
                     df[,colnames(groups)[i]] <- unique(groups[,i])
                     attr(df,"groups") <- c(attr(df, "groups"), colnames(groups)[i])
                   }
                 }
                 df
               })

    # for  each data frame in weightes merger on the scaling factor determinat of
    # curvature (detQ)
    for(wi in 2:length(weights)) {
      weights[[wi]]$detQ <- NULL
      Zgrps <- attr(zScale[[wi]], "groups")
      weights[[wi]] <- merge(weights[[wi]], zScale[[wi]][,c(Zgrps, "detQ")],by.x="index", by.y=Zgrps)
    }

    # make new omega
    omega <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X))
    omegaFull <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X), full=TRUE)
    # make stubs
    Qi0_ <- list(NULL)
    Qi_ <- list(NULL)
    tmpomega <- list(NULL)
    for( oi in 2:length(omega0)) {
      omg0 <- omega0[[oi]]
      omg1 <- 2*omg0  #increment up to allow while loop to run
       
      # keep looking for the mean
      while( max(abs( (omg1 - omg0) / pmax(abs(omg0),1E-5))) > 1E-3) {
        omg1 <- omg0
        tmpomega_ <- c(tmpomega, list(omg0))
        nzi <- ncol(Z[[oi]]) # number of Z columns at olvl
        f <- param.lnl.quad(y, X, oi, Z, ZFull=ZFull, Qi=Qi, QiFull=QiFull,
                            omega, omegaFull=omegaFull, W=weights, k, qp,
                            covariance_constructor, bobyqa=FALSE, verbose=TRUE, acc0=acc, fast=fast,
                            mappedDefault=FALSE)
        for(ici in 1:ncol(omg0)) {
          f0 <- f(par0, top=FALSE, integralMultiplierExponent=0, integralZColumn=ici)
          f1 <- f(par0, top=FALSE, integralMultiplierExponent=1, integralZColumn=ici)
          # f1 is not logged, f0 is
          omg0[,ici] <- as.numeric(f1/f0)
        }
        # make omega
        omega0p <- c(tmpomega, list(omg0))
        while( length(omega0p) < length(omega0)) {
          omega0p[[length(omega0p)+1]] <- omega0[[length(omega0p)+1]] 
        }
        omega <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X))
        omegaFull <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X), full=TRUE)
        if(!adapt) {
          # no need to loop, this is the BLUE
          omg1 <- omg0
        }
      }
      if(verb & adapt) {
        cat("BLUE estimates:\n")
        print(omg0)
      }
      # move the integration points to the new BLUE estiamtes and update
      if(adapt) {
        omg0Full <- buildOmega(omega0=tmpomega_, groups=groups, nrowX=nrow(X), full=TRUE)
        derivatives <- genD(adapterLnL(y, X, levels, Z, ZFull, weights, k, qp,
                                       covariance_constructor, omega,
                                       omg0Full,
                                       tmpomega_, par0, verb, Qi, QiFull, oi, acc),
                            rep(0,sum(unlist(nz)[1:oi], na.rm=TRUE)))
        d2 <- derivatives$D[,-(1:nzi),drop=FALSE] # remove first derivatives
        drv <- d2
        # assumes there are exactly two levels
        Qi0_[[oi]] <- lapply(1:nrow(drv), function(i) { 
          scaleQuadPoints(drv[i,], nzi)
        })
        map <- groups[,oi-1]
        umap <- unique(map)
        Qi_[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
        for(i in 1:nrow(X)) {
          Qi_[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0_[[oi]][[(1:length(umap))[map[i]==umap] ]]
        } 
        QiFull[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
        for(i in 1:nrow(X)) {
          QiFull[[oi]][1:nz,(i-1)*nz+1:nz] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
        }
      }
      # now, actually add this to the tmpomega
      tmpomega <- c(tmpomega, list(omg0))
      omg0Full <- buildOmega(omega0=tmpomega, groups=groups, nrowX=nrow(X), full=TRUE)
    }
    if(adapt) {
      return(list(omega0=tmpomega, omega=omega, omegaFull=omg0Full, Qi0=Qi0_, Qi=Qi_, QiFull=QiFull))
    } else {
      return(tmpomega)
    }
  }
}


# finds the maximum of the likelihood surface and (approximate) SD, per "random effect"
# this function can find these without input estimates. This is used for adapting
MAP <- function(groups, y, X, levels, Z, ZFull, weights, k, qp, covariance_constructor, verbose, nlmevar, nz, acc) {
  ####################################################
  #                        Outline:                  #   
  #        1) Find points for MAP estimation         #
  #        2) MAP estimate                           #
  #         3) Estimate variance                     #
  ####################################################
  
  # when you call adapter, it returns this function that can be called with omega0 and par0
  # and will find a new MAP and SD for each RE
  function(omega0, par0, verb=verbose) {
    #####################################################
    #        1) Find points for MAP estimation         #
    #####################################################
    
    # make omega, which is omega0 at the level one  data length
    omega <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X))
    omegaFull <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X), full=TRUE)
    # f returns the log(likelihood * posterior), per group
    # use Newton's method, per group, to find the maximum of the a posterori surface
    Qi0 <- list(NULL)
    Qi <- list(NULL)
    QiFull <- list(NULL)
    tmpomega <- list(NULL)
    tmpomegaFull <- list(NULL)
    for( oi in 2:length(omega0)) {
      # Over all levels >1 
      omg0 <- omega0[[oi]]
      omg1 <- 1E20*(omg0+1E-15)
      nzi <- nz[[oi]]
      iter <- 1
      while( iter < 25 & max(abs( (omg1 - omg0) / pmax(abs(omg0),1E-5))) > 1E-3) { # iterate while Omega continues to change more than threshold
        iter <- iter + 1
        omg1 <- omg0
        tmpomega_ <- c(tmpomega, list(omg0))
        toF <- buildOmega(omega0=tmpomega_, groups=groups, nrowX=nrow(X), full=TRUE)
        derivatives <- genD(adapterLnL(y, X, levels, Z, ZFull, weights, k, qp,
                                       covariance_constructor, omega, toF,
                                       tmpomega_, par0, verb, Qi, QiFull, oi, acc),
                            rep(0,nz[oi], na.rm=TRUE))
        # this is just second derivatives
        d1 <- derivatives$D[,(1:nzi),drop=FALSE]
        d2 <- derivatives$D[,-(1:nzi),drop=FALSE] # remove first derivatives

        # wrap the second derivatives into matrixes
        d2w <- lapply(1:nrow(d2), function(i) { wrap(d2[i,],nzi)})
        
        ######################################
        #        2) MAP estimate            #
        ######################################
        
        # this is the Newton step, per group
        omg0 <- lapply(1:nrow(d2),
                       function(i) {
                         omg0[i,] - 0.5 * solve(d2w[[i]]) %*% d1[i,]
                       })
        
        # this recudes that to a vector
        omg0 <- t(do.call(cbind, omg0))
        # make omega
        omega0p <- c(tmpomega, list(omg0))
        while( length(omega0p) < length(omega0)) {
          omega0p[[length(omega0p)+1]] <- omega0[[length(omega0p)+1]] 
        }
        omega <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X))
        omegaFull <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X), full=TRUE)
        
      } # end while( max(abs( (omg1 - omg0) / pmax(abs(omg0),1E-5))) > 1E-3)
      if(verb) {
        cat("MAP estimates:\n")
        print(omg0)
      }


      #####################################################
      #         3) Estimate variance (Fishers I)          #
      #####################################################

      drv <- d2
      # assumes there are exactly two levels
      Qi0[[oi]] <- lapply(1:nrow(drv), function(i) {
        ss <- scaleQuadPoints(drv[i,], nzi)
        for(j in 1:nrow(ss)) {
          if(ss[j,j] > abs(omg0[i,j])) {
            ss[j,j] <- sqrt(ss[j,j]^2 + omg0[i,j]^2)
            omg0[i,j] <<- 0
          }
        }
        ss
      })
      map <- groups[,oi-1]
      umap <- unique(map)
      nzi <- ncol(Z[[oi]]) # number of Z columns at olvl
      Qi[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(weights[[oi-1]]))
      for(i in 1:nrow(weights[[oi-1]])) {
        Qi[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
      }
      QiFull[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
      for(i in 1:nrow(X)) {
        QiFull[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
      }
      # now, actually add this to the tmpomega
      tmpomega <- c(tmpomega, list(omg0))
      tmpomegaFull <- omegaFull
      
      if(oi < length(omega0)) {
        # add detQ
        df <- data.frame(detQ=sapply(Qi0[[oi]],det))
        # add group labels
        group_name <- colnames(groups)[oi-1]
        for(i in 1:length(group_name)) {
          df[,group_name[i]] <- unique(groups[,group_name[i]])
        }
        # add detQ to weights for later iterations to use
        weights[[oi]] <- merge(weights[[oi]],df[,c(group_name, "detQ")],by.x="index", by.y=group_name)
      }
    }
    list(omega0=tmpomega, omega=omega, omegaFull=tmpomegaFull, Qi0=Qi0, Qi=Qi, QiFull=QiFull)
  }
}


# helper functions for adapter:

# wrap second derivitive function (genD) output vector into a matrix
wrap <- function(d2, nz) {
  dm <- matrix(0,nrow=nz, ncol=nz)
  for(i in 1:nz) {
    dm[i:nz,i] <- dm[i,i:nz] <- d2[1:(nz-i+1)]
    d2 <- d2[-(1:(nz-i+1))]
  }
  dm
}

# turn the genD output into the Cholesky of the inverse Fisher information
scaleQuadPoints <- function(d2, nz){
  dm <- wrap(d2, nz)
  solved<- solve(-1*dm)
  tryCatch(res <- chol(solved),
           error= function(e) {
             cat("Working on second inverse second derivative matrix\n")
             print(solved)
             stop("non-negative second derivatives with respect to the random effects at their ostensible maxima")
           })
  res
}

# turn omega0 (which has the same number of rows as there are units at each level
# into omega (which has the same number of rows as the X data) 
buildOmega <- function(omega0, groups, nrowX, full=FALSE) {
  omega <- list(NULL)
  oind <- 1
  for(o0i in 2:length(omega0)) {
    omega0i <- as.matrix(omega0[[o0i]])
    res <- matrix(0, nrow=nrowX, ncol=ncol(omega0i))
    noind <- ncol(omega0i) # number of columns to copy over
    map <- groups[,o0i-1]
    umap <- unique(map)
    for(i in 1:length(umap)) {
      # for each unique level
      for(oindi in 1:noind) {
        # copy over data by column (oindi) on omega0i
        # effectively duplicated each value of omega0  a number of times coresponding to how many obs there are in that group
        res[which(map==umap[i]),oindi] <- omega0i[i,oindi]
      }
    }
    if(o0i > 2 & !full) {
      res <- res[!duplicated(groups[,o0i-2]),]
    }
    omega <- c(omega, list(res))
    oind <- oind + noind
  } # closes main for loop for(o0i in 2:length(omega0))
  omega
}

# function that returns the likelihood, by group, evaulated at a point
adapterLnL <- function(y, X, levels, Z, ZFull, weights, k, qp, covariance_constructor, omega, omegaFull, omega0, par0, verbose, Qi, QiFull, olvl, acc) {
  # olvl is the level we are optimizing now
  function(par, long=FALSE) {
    # notice par here is a new displacement in the random effect location
    # not par0, which is a set of mixed model parameters
    yp <- y
    o0 <- omega0
    nzi <- 0 # number of Z columns at olvl
    for(i in 1:olvl) {
      if(!is.null(Z[[i]])) {
        ki <- ncol(Z[[i]])
        if(i == olvl) {
          nzi <- ki
        }
        if(ki >= 1) {
          # keep track of update in omega for calculation of prior prob
          # yp is moved by omega, but also to account for an addition change 
          # controlled by par
          zAdjust <- apply(ZFull[[i]] * omegaFull[[i]],1,sum)
          if(olvl == i) {
            zAdjust <- zAdjust + ZFull[[i]] %*% par[1:ki]
            for(kii in 1:ki) {
              o0[[i]][,kii] <- o0[[i]][,kii] + par[kii]
            }
            par <- par[-(1:ki)]
          }
          # for olvl > 2, zAdjust has the wrong number of rows and need to be mapped
          # back to yp using the weights
          yp <- yp - zAdjust
        } # ends if  if(ki >= 1) 
      } # ends  if(!is.null(Z[[i]]))
    } #ends for(i in 1:olvl)
    # evaluate the likelihood
    beta <- par0[1:k]
    parC <- covariance_constructor(par0[-(1:k)])

    # Set curavtures Q for lower level to matrix of 0 
    # This is becasue the curvature is 0 at the maximum 
    # used for calculated likelihood by groups. 
    Qi_ <- matrix(0, nrow=nzi, ncol=nzi*nrow(weights[[olvl-1]])) # this needs to be the Qi for lower levels
    Qi__ <- c(Qi, list(Qi_))
    QiFull_ <- matrix(0, nrow=nzi, ncol=nzi*nrow(X)) # this needs to be the Qi for lower levels
    QiFull__ <- c(Qi, list(QiFull_))

    #use helper funciton to return the likelihood contribution of teach group 
    loglikelihoodByGroup <- calc.lin.lnl.quad(y=yp, yhat=X %*% beta, level=olvl,
                                              Z, Qi=Qi__,
                                              omega=lapply(omega, function(omegai) {0*omegai}),
                                              W=weights, C=parC, qp, top=FALSE,
                                              atPoint=TRUE, verbose=verbose,
                                              acc=acc, ZFull=ZFull,
                                              omegaFull=omegaFull, QiFull=QiFull__)
    Cl <- parC[[olvl]]
    # apply multivatiate normal pdf to get posterior for each gorup 
    posteriorByGroup <- apply(o0[[olvl]], MARGIN=1, function(p) {
      mvnpdfC(as.matrix(p), rep(0, length = length(p)), varcovM=Cl%*%t(Cl), Log=TRUE)
                        })
    if(long) {
      return(list(res=loglikelihoodByGroup+posteriorByGroup,
                  loglikelihoodByGroup=loglikelihoodByGroup,
                  posteriorByGroup=posteriorByGroup))
    }
    #final return 
    loglikelihoodByGroup + posteriorByGroup
  } #ends  funciton to be return function(par, long=FALSE) 
} 
