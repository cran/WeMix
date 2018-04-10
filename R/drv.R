# calculate a numerical central derivative 
d <- function(f, par, delta=1E-5) {
  res <- rep(NA, length(par))
  for(i in 1:length(par)) {
    d <- rep(0,length(par))
    d[i] <- delta/2
    fp <- f(par + d)
    fm <- f(par - d)
    res[i] <- (fp - fm) / delta
  }
  res
}

#' @importFrom numDeriv hessian
d2<-function(lnlfn, par){
  hessian(lnlfn, par, method="Richardson")
}

# This function uses numDeriv package to calculate both grad and hessian using genD
getHessianAndGrad <- function(func, x){
  # builds on math from numDeriv 
  # see citation in vignette to Gilbert & Ravi.(2016)
  
  args <- list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), 
               r=4, v=2, show.details=FALSE) # default
  
  D <- genD(func, x, method="Richardson", method.args=args)$D
  if(1!=nrow(D)) stop("Derivative from getHessianAndGrad function has too many rows")
  H <- diag(NA,length(x))
  u <- length(x)
  for(i in 1:length(x))
    for(j in 1:i){ 
      u <- u + 1
      H[i,j] <- D[,u]
    }
  
  H <- H + t(H)
  diag(H) <- diag(H)/2
  H
  return(list(grad = D[1:length(x)],
              hessian = H))
}
