r_irlba <- function(X, nu, work=nu+7, maxit=1000, tol=1e-5, eps=1e-10, svtol=tol) {
  stopifnot(work>nu)
  IRLB(X, nu, work, maxit, tol, eps, svtol)
}

#' @import irlba
r_orthog <- function(x, y) {
  if (missing(y))
    y <- runif(nrow(x))
  y <- matrix(y)
  xm <- nrow(x)
  xn <- ncol(x)
  yn <- ncol(y)
  stopifnot(nrow(y)==xm)
  stopifnot(yn==1)
  initT <- matrix(0, xn+1, yn+1)
  ORTHOG(x, y, initT, xm, xn, yn)
}
