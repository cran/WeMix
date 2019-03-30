#discrepency function for calculating least squares solution. Follows from Bates et al. 2015 eq 14 and 15
discf <- function(y, Zt, X, lambda, u, Psi12, W12, b) { # return ehatT * ehat
  if(any(is.na(b))) {
    X <- X[,!is.na(b),drop=FALSE]
    b <- b[!is.na(b)]
  }
  wres <- W12 %*% (y - Matrix::t(Zt) %*% lambda %*% u - X %*% b) # residual
  ures <- Psi12 %*% u # augmented ehat
  as.numeric(sum(wres^2) + sum(ures^2))
}

# get qr, first making sure to use sparse methods so we can use sparse functions and avoid many if/else statements
qr_0 <- function(X, allowPermute=FALSE) { 
  # Matrix::qr does not do well in this case, use base::qr
  if(nrow(X) < ncol(X)) {
    X <- as(X, "matrix")
    return(base::qr(X))
  }
  # try to use Matrix
  qr1 <- Matrix::qr(X)
  if(allowPermute || inherits(qr1,"qr")) {
    # if allowPermute (allow Matrix::qr to use permutations) or this is a base::qr, just return
    return(qr1)
  }
  # do not allow permutations. Test if they happened. If they did, use base::qr instead
  if( any(range( order(qr1@q) - 0:max(qr1@q)) > 0) ) {
    X <- as(X, "matrix")
    return(base::qr(X))
  }
  # base case, return the non-permuted Matrix::qr result
  return(qr1)
}

qr_s <- function(X, allowPermute=FALSE) { 
  X <- as(X, "matrix")
  return(base::qr(X))
}

# directly return R from a QR
qr_qrr <- function(X) { 
  if(nrow(X) < ncol(X)) {
    X <- as(X, "matrix")
    return(base::qr.R(base::qr(X)))
  }
  qr1 <- Matrix::qr(X)
  if(inherits(qr1,"qr")) {
    return(base::qr.R(qr1))
  }
  # check for permutations and use base if there are some
  if( any(range( order(qr1@q) - 0:max(qr1@q)) > 0) ) {
    return(base::qr.R(base::qr(X)))
  }
  # we can use Matrix:qr
  return(Matrix::qrR(qr1, backPermute=FALSE))
}


qrr_ <- function(QR) {
  if(inherits(QR, "qr")) {
    return(base::qr.R(QR))
  } else {
    return(Matrix::qrR(QR, backPermute=FALSE))
  }
}

# directly return R from a chol of AtA
getchr <- function(A) {
  suppressWarnings(m0 <- as(Matrix::chol(Matrix::t(A) %*% A, pivot=FALSE), "Matrix"))
  colnames(m0) <- colnames(A) # allowable because pivot=FALSE
  return(m0)
}

chr_ <- function(A, qr) {
  m <- tryCatch(getchr(A),
                error=function(e) {
                  if(inherits(qr,"qr")) {
                    return(base::qr.R(qr))
                  } else {
                    A <- as(A, "matrix")
                    return(base::qr.R(base::qr(A)))
                  }
                })
  return(m)
}

rkfn <- function(qr_R) {
  d <- Matrix::diag(qr_R)
  tol <- length(d) * .Machine$double.eps
  sum(abs(d) > max(abs(d)) * tol)
}


#' The main function which calculates the analytic solution to the linear mixed effects model. 
#' @param weights level-1 weights
#' @param y outcome measure. 
#' @param X the X matrix.
#' @param Zlist, a list of matrixes with the Z values for each level. 
#' @param Zlevels, the level corresponding to each matrix in Zlist. 
#' @param weights a list of unconditional weights for each model level. 
#' @param weightsC a list of conditional weights for each model level. 
#' @param groupID a matrix containing the group ids for each level in the model. 
#' @param lmeVardf a dataframe containing the variances and covariance of the random effects, in the same format as returned from lme. 
#' @importFrom Matrix Diagonal Matrix bdiag
#' @importFrom methods as .hasSlot
#' @keywords internal
analyticSolve <- function(y, X, Zlist, Zlevels, weights, weightsC=weights, groupID, lmeVarDF, analyticJacobian=FALSE, forcebaseQR=FALSE, qr_=qr_0, v0) {
  if(!inherits(groupID, "matrix")) {
    stop(paste0("Variable", dQuote(groupID)," must be a matrix with IDs at a level per column."))
  }
  if(!inherits(y, "numeric")) {
    stop(paste0(dQuote("y"), " must be a numeric vector."))
  }
  ny <- length(y)
  X <- as.matrix(X)
  nX <- nrow(X)
  if(nX != ny) {
    stop(paste0("Length of the ", dQuote("y"), " vector and the ", dQuote("X"), " matrix must agree."))
  }
  if(length(Zlist) != length(Zlevels)) {
    stop(paste0("Length of the ", dQuote("Zlist"), " and ", dQuote("Zlevels"), " must agree."))
  }
  nW1 <- length(weights[[1]])
  if(nW1 != nX) {
    stop(paste0("Number of rows in ", dQuote("weights[[1]]") , " must agree with number of rows in ", dQuote("X"),"."))
  }
  nWc1 <- length(weightsC[[1]])
  if(nWc1 != nX) {
    stop(paste0("Number of rows in ", dQuote("weightsC[[1]]") , " must agree with number of rows in ", dQuote("X"),"."))
  }
  ngid <- nrow(groupID)
  if(ngid != nX) {
    stop(paste0("Number of rows in ", dQuote("groupID") , " must agree with number of rows in ", dQuote("X"),"."))
  }
  
  # get the number of groups per level
  nc <- apply(groupID, 2, function(x) {length(unique(x))} )
  
  # for level 1 cast Z for lowest level as a matrix and transpose to get Zt
  Zt <- Matrix::t(as(Zlist[[1]], "Matrix"))
  
  # for level, build PsiVec, the diagonal of the Psi  (weights) matrix
  PsiVec <- rep(weights[[Zlevels[1]]], ncol(Zlist[[1]])/nc[1])
  
  # Assemble the Zt matrix and psi Vector for levels >1 
  ZColLevels <- rep(Zlevels[1], ncol(Zlist[[1]]))
  if(length(Zlist) > 1) {
    for(zi in 2:length(Zlist)) {
      Zt <- rbind(Zt, Matrix::t(as(Zlist[[zi]], "Matrix")))
      ZColLevels <- c(ZColLevels, rep(Zlevels[zi], ncol(Zlist[[zi]])))
      PsiVec <- c(PsiVec, rep(weights[[Zlevels[zi]]], ncol(Zlist[[zi]])/nc[Zlevels[zi]-1]))
    }
  }
 
  # get unit weights
  W0 <- weights[[1]]
  Psi <- Diagonal(n=nrow(Zt), x=PsiVec) #matrix of weights at level one
  Psi12 <- Diagonal(n=nrow(Zt), x=sqrt(PsiVec)) #matrix of square root of level one weights
  W <- Diagonal(x=(W0))
  W12 <- Diagonal(x=sqrt(W0))
  
  # calculate  conditional weights matrix where level one weights are scaled by level two wieghts
  W12cDiag <- W0
  L1groups <- unique(groupID[,1])
  for(gi in 1:length(L1groups)) {
    rowsGi <- groupID[,1] == L1groups[gi]
    W12cDiag[rowsGi] <- sqrt(W12cDiag[rowsGi] / weights[[2]][gi])
  }
  W12c <- Diagonal(x=W12cDiag)
  Z <- Matrix::t(Zt)
  options(Matrix.quiet.qr.R=TRUE) #supress extra print outs 

  M0 <- Matrix(data=0, ncol=ncol(X), nrow=nrow(Zt))


  # build pseudo-Delta (not actual theta)
  iDelta <- Delta <- list()
  levels <- max(lmeVarDF$level)
  for (l in 2:levels){
    #set up matrix of zeros with rows and columns for each random effect
    #number of random effect is the number of variance terms at that level (not including covariance terms )
    n_ref_lev  <- nrow(lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),] ) 
    lambda_i  <- matrix(0,nrow=n_ref_lev,ncol=n_ref_lev)
    row.names(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),"var1"]
    colnames(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),"var1"]
    
    #get only v for this level
    group <- lmeVarDF[lmeVarDF$level==l ,"grp"][1]
    v_lev <- v0[grep(paste0("^",group,"."),names(v0))]
    
    
    #fill in lambda_i from theta using names
    for (vi in 1:length(v_lev)){
      row_index <- strsplit(names(v_lev[vi]),".",fixed=TRUE)[[1]][2]
      col_index <- ifelse(length(strsplit(names(v_lev[vi]),".",fixed=TRUE)[[1]])>2,strsplit(names(v_lev[vi]),".",fixed=TRUE)[[1]][3],row_index)
      lambda_i[row_index,col_index] <- v_lev[vi]
    }
    diag(lambda_i)[diag(lambda_i)==0] <- .Machine$double.eps #set any 0s to smallest non zero number to enable solve in next step
    iDelta[[l]] <- lambda_i
    Delta[[l]] <- solve(lambda_i)
  }
  nozero <- function(x) { any(x != 0) }

  ZiAl <- list()
  for(level in 1:ncol(groupID)) {
    groups <- unique(groupID[,level])
    Deltai <- Delta[[level+1]]
    topLevel <- level == ncol(groupID)
    Zl <- Z[,ZColLevels==(level+1),drop=FALSE] # level l Z
    for(gi in 1:length(groups)) {
      # this works for 3 level, will need to consider selection by row for >3 
      Zrows <- groupID[,level]==groups[gi]
      # filter to the rows for this group, also filter to just columns
      # for this level of the model
      Zi <- Zl[Zrows,,drop=FALSE] # Zi associated with the random effect (~ Z_g in the specs)
      # within columns for this level, only the non-zero columns
      # regard this unit, so filter it thusly
      # the below code uses the pointers from the sparse matrix to identify non zero rows
      # it is equivalent to
      # Zcols <- apply(Zi,2,nozero)
  
      if(.hasSlot(Zi,"p")){
        z_pointers <- Zi@p
        len_pointers <- length(z_pointers)
        Zcols <- which(z_pointers[1:(len_pointers-1)]  - z_pointers[2:len_pointers] !=  0 ) 
      } else {
        Zcols <- apply(Zi,2,nozero)
      }

      Zi <- Zi[,Zcols,drop=FALSE]

      if(level == 1 || level < ncol(groupID)) {
        if(!topLevel) {
          lp1 <- unique(groupID[Zrows,level+1])
          if(length(lp1) > 1) {
            stop("Not a nested model.")
          }
          qp1 <- nrow(Delta[[level+2]])
        }
        qi <- nrow(Deltai) 
        Zi <- Z[Zrows,,drop=FALSE]
        # within columns for all levels, only the non-zero columns
        # regard this unit, filtered using the pointers from sparse matrix
        # if they are available. this is equivalent to
        # Zcols <- apply(Zi,2,nozero)
        if(.hasSlot(Zi,"p")){
          z_pointers <- Zi@p
          len_pointers <- length(z_pointers)
          Zcols <- which(z_pointers[1:(len_pointers-1)]  - z_pointers[2:len_pointers] !=  0 ) 
        } else {
          Zcols <- apply(Zi,2,nozero)
        }
        ZiA <- W12cDiag[Zrows] * Zi[,Zcols,drop=FALSE]
        ZiAl <- c(ZiAl, list(ZiA))
      }
    }
  }
  
  function(v, verbose=FALSE, beta=NULL, sigma=NULL, robustSE=FALSE, returnBetaf=FALSE, getGrad=FALSE) {
    beta0 <- beta
    sigma0 <- sigma
    omegaVec <- v

    #Create list delta and lambda Matrixes for each level
    # iDelta is the inverse Delta (lambda pre kronecker) for forming VC estimates in the end
    iDelta <- Delta <- list()
    lambda_by_level <- list()
    
    levels <- max(lmeVarDF$level)
   
    for (l in 2:levels){
      #set up matrix of zeros with rows and columns for each random effect
      #number of random effect is the number of variance terms at that level (not including covariance terms )
      n_ref_lev  <- nrow(lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),] ) 
      lambda_i  <- matrix(0,nrow=n_ref_lev,ncol=n_ref_lev)
      row.names(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),"var1"]
      colnames(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),"var1"]
      
      #get only v for this level
      group <- lmeVarDF[lmeVarDF$level==l ,"grp"][1]
      v_lev <- v[grep(paste0("^",group,"."),names(v))]
      
      
      #fill in lambda_i from theta using names
      for (vi in 1:length(v_lev)){
        row_index <- strsplit(names(v_lev[vi]),".",fixed=TRUE)[[1]][2]
        col_index <- ifelse(length(strsplit(names(v_lev[vi]),".",fixed=TRUE)[[1]])>2,strsplit(names(v_lev[vi]),".",fixed=TRUE)[[1]][3],row_index )
        lambda_i[row_index,col_index] <- v_lev[vi]
      }
      diag(lambda_i)[diag(lambda_i)==0] <- .Machine$double.eps #set any 0s to smallest non zero number to enable solve in next step
      iDelta[[l]] <- lambda_i
      Delta[[l]] <- solve(lambda_i)
      
      #assemble blockwise lambad by level matrix 
      lambda_by_level[[l]] <- kronecker(lambda_i,diag(1,unique(lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),"ngrp"])))
    }
  
    #Create lambda matrix from diagonal of all lamdas 
    lambda <- bdiag(lambda_by_level[2:length(lambda_by_level)]) # vignette equation 21

    yaugmented <- c(as.numeric(W12 %*% y), rep(0, nrow(Zt))) # first matrix in vingette equation 25 and 60 
    
    A <- rbind(cbind(W12 %*% Z %*% lambda, W12 %*% X),
               cbind(Psi12,M0)) # second matrix in equation (60)/(25). 
    qr_a <- qr_(A, allowPermute=TRUE)
    suppressWarnings(ub <- Matrix::qr.coef(qr_a, yaugmented)) # equation (26)/ stack vector (u: top part for fixed effects coef estimate, b: bottom part for random effects vector (i.e. of length n_group * n_random_effects)

    # get the gradients
    if(getGrad) {
      b <- ub[-(1:ncol(Z))]
      u <- ub[1:ncol(Z)]
      R <- chr_(A, qr_a)
      Rrows <- nrow(R)-length(b):1+1
      Rrows <- Rrows[Rrows > length(u)]
      RX <- R[Rrows,ncol(R)-length(b):1+1,drop=FALSE]
      db <- beta0 - b
      db <- db[colnames(RX)]
      db[is.na(db)] <- beta0[is.na(db)]
      res0 <- res <- 0 * db
      for(bi in 1:length(b)) {
        ei <- 0 * res
        ei[bi] <- 1
        #(RX*RX) is elementwise multiplication
        res[bi] <- -1*(db[bi]) * sum((RX*RX) %*% ei)/sigma^2
      }
      return(res)
    }

    # get ln of determinant of Lz matrix (Bates et al. 2015 notation) / calculation of alpha
    # or the final product in Pinheiro and Bates 2.13 (formula used for weights)
    lndetLz <- 0 # lndetLz is alpha in the specs
    lndetLzg <- list() # by top level group lnldetLz
    # used to find columns with no-non-zero elements
    nozero <- function(x) { any(x != 0) }
    for(level in 1:ncol(groupID)) {
      groups <- unique(groupID[,level])
      Deltai <- Delta[[level+1]]
      sDeltai0 <- as(Delta[[level+1]], "sparseMatrix")
      # R11 notation from Bates and Pinheiro, 1998
      R11 <- list()
      topLevel <- level == ncol(groupID)
      for(gi in 1:length(groups)) {
        if(level == 1 || level < ncol(groupID)) {
          # this may work for 3 level, will need to consider selection by row for >3 
          Zrows <- groupID[,level]==groups[gi]
          # get Zi
          ZiA <- ZiAl[[gi]]
          if(!topLevel) {
            lp1 <- unique(groupID[Zrows,level+1])
            qp1 <- nrow(Delta[[level+2]])
          }
          qi <- nrow(Deltai) 
          #shape delta, used to rbind like this:
          #Delta0 <- Matrix(data=0, nrow=qi, ncol=ncol(ZiA)-ncol(Deltai))
          #Delta0 <- sparseMatrix(i=integer(0), j=integer(0), dims=c(qi, ncol(ZiA)-ncol(Deltai)))
          #ZiA <- rbind(ZiA, cbind(Deltai, Delta0))
          # expand sDeltai as if we rbinded, much faster
          sDeltai <- sDeltai0          
          sDeltai@Dim <- c(sDeltai@Dim[1], sDeltai@Dim[2] + ncol(ZiA)-ncol(Deltai))
          sDeltai@Dimnames[[2]] <- rep("", sDeltai@Dim[2])
          sDeltai@p <- c(sDeltai@p, rep(max(sDeltai@p),ncol(ZiA)-ncol(Deltai)))

          ZiA <- rbind(ZiA, sDeltai)
          #in this section we decompose Zi * A using qr decomposition, following equation 43 in the vingette.
          ZiAR <- qr_qrr(ZiA) # this is faster than getchr(ZiA)
          R22 <- ZiAR[1:qi, 1:qi, drop=FALSE]
          lndetLzi <- weights[[level+1]][gi]*(- as.numeric(determinant(Deltai)$mod) + sum(log(abs(Matrix::diag(R22)[1:ncol(Deltai)])))) # equation (92)
          # save R11 for next level up, if there is a next level up
          if(!topLevel) {
            R11i <- sqrt(weightsC[[level+1]][gi])*ZiAR[-(1:qi),qi+1:qp1, drop=FALSE]
            # load in R11 to the correct level
            if(length(R11) >= lp1) {
              R11[[lp1]] <- rbind(R11[[lp1]], R11i)
              lndetLzg[[lp1]] <- lndetLzg[[lp1]] + lndetLzi
            } else {
              # there was no R11, so start it off with this R11
              R11[[lp1]] <- R11i
              lndetLzg[[lp1]] <- lndetLzi
            }
          } else {
            lndetLzg[[gi]] <- lndetLzi
          }
        }
        if(level>=2) {
          ZiA <- rbind(pR11[[groups[gi]]], Deltai)
          R <- qr_qrr(ZiA) # probably faster than getchr(ZiA)
          # weight the results
          lndetLzi <- weights[[level+1]][gi]*(- (as.numeric(determinant(Deltai)$mod)) + sum(log(abs(Matrix::diag(R)[1:ncol(Deltai)])))) # equation 92
          lndetLzg[[gi]] <- lndetLzg[[gi]] + lndetLzi
        }
        # update W^{1/2}, for this purpose, it's the conditional weights
        # W12i <- W12c[Zrows,Zrows]
        lndetLz <- lndetLz + lndetLzi
        if(verbose) {
          cat("gi=", gi, " lndetLzi=", lndetLzi, " lndetLz=", lndetLz,"w=", weights[[level+1]][gi], "\n")
        }
      }
      if(level < ncol(groupID)) {
        pR11 <- R11
      }
    }
    # this is the beta vector
    b <- ub[-(1:ncol(Z))]
    if(!is.null(beta0)) {
      if(length(b) != length(beta0)) {
        stop(paste0("The argument ", dQuote("beta"), " must be a vector of length ", length(b), "."))
      }
    }
    # and the random effect vector
    u <- ub[1:ncol(Z)]
  
    
    # wrap the random effect vector to matrix format so that each row regards one group
    bb <- list()
    vc <- matrix(nrow=1+length(v), ncol=1+length(v))
    u0 <- 0
    for(li in 2:length(lambda_by_level)) {
      lli <- lambda_by_level[[li]]
      ni <- nrow(lli)
      bb[[li]] <- lli %*% u[u0 + 1:ni]
      u0 <- u0 + ni
      vc[li,li] <- 1/(ni) * sum((bb[[li]])^2) # mean already 0
    }

    # the discrepancy ||W12(y-Xb-Zu)||^2 + ||Psi(u)||^2
    discrep <- discf(y, Zt, X, lambda, u, Psi12, W12, b)
    # the R22 matrix, bottom right of the big R, conforms with b
    R <- chr_(A, qr_a)
    Rrows <- nrow(R)-length(b):1+1
    Rrows <- Rrows[Rrows > length(u)]
    RX <- R[Rrows,ncol(R)-length(b):1+1,drop=FALSE]
    nx <- sum(W0)
    # residual
    sigma <- sqrt(discrep / nx)
    dev <- 0
    if(!is.null(sigma0)) {
      sigma <- sigma0
    }
    dev <- 2*as.numeric(lndetLz) + nx*log(2*pi*(sigma^2)) + discrep/sigma^2
    # add R22 term if beta is not beta-hat
    if(returnBetaf) {
      f <- function(beta0) {
        db <- beta0 - b
        # if there is an NA, then use the beta value--that is, assume betahat is 0
        db[is.na(db)] <- beta0[is.na(db)]
        # RX will be reordered if X has all 0 columns, db needs to be rearanged similarly
        db <- db[colnames(RX)]
        # Matrix::qr.R makes bad colnames, fix that
        db[is.na(db)] <- 0
        # deviance
        dev <- dev + sum( (RX %*% db)^2 )/sigma^2
        res <- return(dev/-2)
      }
      return(f)
    }
    # add in term for beta not beta-hat
    if(!is.null(beta0)) {
      # difference in beta and betahat (b)
      db <- beta0 - b
      # if there is an NA, then use the beta value--that is, assume betahat is 0
      db[is.na(db)] <- beta0[is.na(db)]
      # RX will be reordered if X has all 0 columns, db needs to be rearanged similarly
      db <- db[colnames(RX)]
      # deviance
      dev <- dev + sum( (RX %*% db)^2 )/sigma^2
    }
    #Variance of Beta, from Bates et.al. 2015
    # qr matrix rank, tolerance value and method based on Matrix::rankMatrix API
    if(rkfn(R) == ncol(A)) {
      RXi <- Matrix::solve(RX)
      cov_mat <- (sigma^2) * RXi %*% Matrix::t(RXi)
      varBeta <- Matrix::diag(cov_mat)
    } else {
      cov_mat <- NULL
      varBeta <- rep(NA, length(b))
    }

    sigma <- ifelse(is.null(sigma0), sqrt(discrep / nx), sigma0)

    res <- list(b=b, u=u, ranef=bb, discrep=discrep, sigma=sigma, dev=dev,
                lnl=dev/-2, theta=v, varBeta=varBeta, vars=NULL, cov_mat=cov_mat,
                lndetLz=lndetLz, iDelta=iDelta)
    if(robustSE) {
      #Calculate robust standardize effect
      # based on the sandwich estimator of Rabe-Hesketh and Skrondal 2006
      # this whole calculation focuses on calculating the liklihood at the top group level 
      fgroupID <- groupID[,ncol(groupID)] 
      uf <- unique(fgroupID) # unique final groupIDs
      # store lnl to check later
      lnli2<- lnli <- vector(length=length(uf))
      Jacobian <- matrix(NA, nrow=length(uf), ncol=length(b)+length(v)+1)
      bwiL <- list()
      #this code block seperates wieghts into the set of weights belonging to each top level group 
      # components of discrep
      wres <- W12 %*% (y - Matrix::t(Zt) %*% lambda %*% u - X %*% b) # residual
      ures <- Psi12 %*% u # augmented ehat
      for(gi in 1:length(uf)) {
        sgi <- fgroupID == gi # subset for group gi
        weightsPrime <- list(weights[[1]][sgi])
        weightsCPrime <- list(weightsC[[1]][sgi])
        for(i in 2:length(weights)) {
          theseGroups <- unique(groupID[sgi,i-1])
          weightsPrime[[i]] <- weights[[i]][theseGroups]
          weightsCPrime[[i]] <- weightsC[[i]][theseGroups]
        }
        condenseZ <- function(z) {
          res <- z[sgi,,drop=FALSE]
          res[,apply(abs(res),2,sum)>0, drop=FALSE]
        }
        groupIDi <- groupID[sgi,,drop=FALSE]
        lmeVarDFi <- lmeVarDF
        lmeVarDFi$ngrp[lmeVarDFi$level==1] <- sum(sgi)
        for(i in 2:length(weights)) {
          lmeVarDFi$ngrp[lmeVarDFi$level==i] <- length(unique(groupIDi[,i-1]))
        }
        #calculate the group level likelihood by applying the function to only the x or y within the group indexed by sgi
lnli2[gi] <- sum(wres[sgi]^2)/sigma^2/-2 + (ures[gi]^2)/sigma^2/-2

        bwi <- analyticSolve(y=y[sgi], X[sgi,,drop=FALSE],
                  Zlist=lapply(Zlist, condenseZ),
                  Zlevels=Zlevels,
                  weights=weightsPrime,
                  weightsC=weightsCPrime,
                  groupID=groupIDi,
                  lmeVarDF=lmeVarDFi,
                  v0=v) 
        tryCatch(lnli[gi] <- bwi(v=v, verbose=verbose, beta=b, sigma=sigma, robustSE=FALSE)$lnl,
                 error= function(e) {
                   lnli[gi] <<- NA
                 })

        # sometimes Matrix::qr tries to solve a singular system and fails
        # normally base::qr works in these cases, so use that instead
        if( abs(lnli2[gi] - lnli[gi]) > 0.1) {
          bwi <- analyticSolve(y=y[sgi], X[sgi,,drop=FALSE],
                               Zlist=lapply(Zlist, condenseZ),
                               Zlevels=Zlevels,
                               weights=weightsPrime,
                               weightsC=weightsCPrime,
                               groupID=groupIDi,
                               lmeVarDF=lmeVarDFi,
                               qr_=qr_s,
                               v0=v) 
           bwiL <- c(bwiL, list(bwi))
           tryCatch(lnli[gi] <- bwi(v=v, verbose=verbose, beta=b, sigma=sigma, robustSE=FALSE)$lnl,
                    error= function(e) {
                      lnli[gi] <<- NA
                    })
        } else {
           bwiL <- c(bwiL, list(bwi))
        }

        # if the function can be evaluated for this group.
        if(!is.na(lnli[gi])) {
          bwiW <- function(bwi, v, sigma, b, ind) {
            bwip <- bwi(v=v, verbose=FALSE,
                        b=b, sigma=sigma,
                        robustSE=FALSE, returnBetaf=TRUE)
            function(par) {
              b[ind] <- par
              bwip(b)
            }
          }
          bwiT <- function(bwi, v, sigma, t, ind) {
            function(par) {
              if(ind > length(v)) {
                sigma <- par
              } else {
                v[ind] <- par
              }
              bwi(v=v, verbose=FALSE,
                  b=b, sigma=sigma, robustSE=FALSE)$lnl
            }
          }
          if(analyticJacobian) {
            Jacobian[gi,1:length(b)] <- bwi(v=v, sigma=sigma, b=b, robustSE=FALSE, getGrad=TRUE)
          } else {
            for(j in 1:length(b)){
              Jacobian[gi,j] <- d(bwiW(bwi, v=v, sigma=sigma, b=b, ind=j), par=b[j])
            }
          }
          for(j in 1:(1+length(v))) {
            v0 <- c(v, sigma)[j]
            Jacobian[gi,length(b)+j] <- d(bwiT(bwi, v=v, sigma=sigma, b=b, ind=j), par=v0)
          }
        } else { # end if(!is.na(lnli[gi]))
          warning("Robust SE downward biased. Some top level units not large enough. Try combining.")
        }
      }
      lnl1 <- dev/-2
      lnliT <- sum(lnli)
      if(!is.na(lnliT)) {
        if( abs(lnl1 - lnliT) > 0.1) {
          # this would be a large difference
          warning("Likelihood estimated at the top group level and summed disagrees with overall likelihood. Standard errrors may not be accurate.")
        }
      }
      J <- matrix(0,ncol=ncol(Jacobian), nrow=ncol(Jacobian))
      nr <- 0
      # a cross product that removes rows (groups) where lnl cannot be evaluated
      for (i in 1:nrow(Jacobian)){
        if(!is.na(lnli[i])) {
          J  <- J  + Jacobian[i,]  %*% t(Jacobian[i,])
          nr <- nr + 1
        }
      }
      J <- (nr/(nr-1))*J
      # just beta part of J
      varBetaRobust <- as(cov_mat %*% J[1:length(b), 1:length(b)] %*% cov_mat , "matrix")
      colnames(J) <- rownames(J) <- c(names(b), names(v), "sigma")
      res <- c(res, list(varBetaRobust=varBetaRobust,
                         seBetaRobust=sqrt(diag(varBetaRobust)), iDelta=iDelta,
                         Jacobian=J))
    }
    return(res)
  }
}

#' a wrapper for analyticSolve that allows it to be called from an optimizer. Takes the same arguments as analyticSolve. 
#' @param weights level-1 weights
#' @param y outcome measure. 
#' @param X the X matrix.
#' @param Zlist, a list of matrixes with the Z values for each level. 
#' @param Zlevels, the level corresponding to each matrix in Zlist. 
#' @param weights a list of unconditional weights for each model level. 
#' @param weightsC a list of conditional weights for each model level. 
#' @param  groupID a matrix containing the group ids for each level in the model. 
#' @param  lmeVardf a dataframe containing the variances and covariance of the random effects, in the same format as returned from lme. 
#' @importFrom Matrix Diagonal Matrix
#' @importFrom methods as
#' @keywords internal
devG <- function(y, X, Zlist, Zlevels, weights, weightsC=weights, groupID, lmeVarDF, v0) {
  bs <- analyticSolve(y=y, X=X, Zlist=Zlist, Zlevels=Zlevels, weights=weights, weightsC=weightsC, lmeVarDF=lmeVarDF,
           groupID=groupID, v0=v0)
  function(v, getBS=FALSE) {
    if(getBS) {
      return(bs)
    }
    bhat <- bs(v)
    return(bhat$dev)
  }
}
