#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define orthog_sig double*, double*, double*, int, int, int

// [[Rcpp::export]]
arma::mat ORTHOG(arma::mat& X, arma::mat& Y, arma::mat& T, int xm, int xn, int yn) {

  static SEXP(*c_orthog)(orthog_sig) = (SEXP(*)(orthog_sig)) R_GetCCallable("irlba", "orthog");
  // Copy preserve R's data
  arma::mat Ycopy = arma::mat(Y.memptr(), Y.n_rows, Y.n_cols);
  c_orthog(X.memptr(), Ycopy.memptr(), T.memptr(), xm, xn, yn);
  return Ycopy;
}


#define irlb_sig double*, void*, int, int, int, int, int, int, int, double, double*, double*, double*, double*, double*, double*, int*, int*, double, int, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double, double*

/* irlb C++ implementation wrapper for Armadillo
* X double precision input matrix
* NU integer number of singular values/vectors to compute must be > 3
* INIT double precision starting vector length(INIT) must equal ncol(X)
* WORK integer working subspace dimension must be > NU
* MAXIT integer maximum number of iterations
* TOL double tolerance
* EPS double invariant subspace detection tolerance
* MULT integer 0 X is a dense matrix (dgemm), 1 sparse (cholmod)
* RESTART integer 0 no or > 0 indicates restart of dimension n
* RV, RW, RS optional restart V W and S values of dimension RESTART
*    (only used when RESTART > 0)
* SCALE either NULL (no scaling) or a vector of length ncol(X)
* SHIFT either NULL (no shift) or a single double-precision number
* CENTER either NULL (no centering) or a vector of length ncol(X)
* SVTOL double tolerance max allowed per cent change in each estimated singular value */
// [[Rcpp::export]]
List IRLB(arma::mat& X,
                 int nu,
                 int work,
                 int maxit=1000,
                 double tol=1e-5,
                 double eps=1e-9,
                 double svtol=1e-5)
{

  int m = X.n_rows;
  int n = X.n_cols;
  int iter, mprod;
  int lwork = 7 * work * (1 + work);


  arma::vec s = arma::randn<arma::vec>(nu);
  arma::mat U = arma::randn<arma::mat>(m, work);
  arma::mat V = arma::randn<arma::mat>(n, work);

  arma::mat V1 = arma::zeros<arma::mat>(n, work); // n x work
  arma::mat U1 = arma::zeros<arma::mat>(m, work); // m x work
  arma::mat  W = arma::zeros<arma::mat>(m, work);  // m x work  input when restart > 0
  arma::vec F  = arma::zeros<arma::vec>(n);     // n
  arma::mat B  = arma::zeros<arma::mat>(work, work);  // work x work  input when restart > 0
  arma::mat BU = arma::zeros<arma::mat>(work, work);  // work x work
  arma::mat BV = arma::mat(work, work);  // work x work
  arma::vec BS = arma::zeros<arma::vec>(work);  // work
  arma::vec BW = arma::zeros<arma::vec>(lwork); // lwork
  arma::vec res = arma::zeros<arma::vec>(work); // work
  arma::vec T = arma::zeros<arma::vec>(lwork);  // lwork
  arma::vec svratio = arma::zeros<arma::vec>(work); // work

 static SEXP(*c_irlb)(irlb_sig) = (SEXP(*)(irlb_sig)) R_GetCCallable("irlba", "irlb");


  c_irlb (X.memptr(), NULL, 0, m, n, nu, work, maxit, 0,
          tol, NULL, NULL, NULL,
          s.memptr(), U.memptr(), V.memptr(), &iter, &mprod,
          eps, lwork, V1.memptr(), U1.memptr(), W.memptr(),
          F.memptr(), B.memptr(), BU.memptr(), BV.memptr(),
          BS.memptr(), BW.memptr(), res.memptr(), T.memptr(),
          svtol, svratio.memptr());
  return List::create(Rcpp::Named("d")=s,
                            Rcpp::Named("u")=U.cols(0, nu-1),
                            Rcpp::Named("v")=V.cols(0,nu-1),
                            Rcpp::Named("iter")=iter,
                            Rcpp::Named("mprod")=mprod);
                            // Rcpp::Named("converged")=conv);
}
