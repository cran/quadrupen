/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_UTILS_H
#define _quadrupen_UTILS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

inline vec  get_lambda1(SEXP LAMBDA1, SEXP N_LAMBDA, SEXP MIN_RATIO, double lmax) {
  vec lambda1 ;
  if (LAMBDA1 != R_NilValue) {
    lambda1  = as<vec>(LAMBDA1)  ;
  } else {
    uword n_lambda = as<uword>(N_LAMBDA) ;
    double min_ratio = as<double>(MIN_RATIO);
    lambda1 = exp10(linspace(log10(lmax), log10(min_ratio*lmax), n_lambda)) ;
  }
  return(lambda1);
}

inline sp_mat get_struct(SEXP STRUCT, double lambda2, vec penscale) {
  uword p = penscale.n_elem;
  sp_mat S;
  if (STRUCT == R_NilValue | lambda2 == 0) {
    S = speye(p,p);
  } else {
    S = as<sp_mat> (STRUCT) ;
  }
  if (lambda2 > 0) {
    // renormalize the l2 structuring matrix according to the l1
    // penscale values, so as it does not interfer with the l2 penalty.
    S = diagmat(sqrt(lambda2)*pow(penscale,-1/2)) * S * diagmat(sqrt(lambda2)*pow(penscale,-1/2)) ;
  }
  return(S) ;
}

vec  cg(mat A, vec b, vec x, double tol) ;
vec pcg(mat A, mat P, vec b, vec x, double tol) ;

void cholupdate(mat &R, mat& XtX) ;

void choldowndate(mat &R, int j) ;

void bound_to_optimal(vec &betaA, mat &xAtxA, vec &xty, vec &grd, double &lambda1, double &lambda2, double &normy, uvec &A, int &monitor, vec &J_hat, vec &D_hat) ;

void add_var_enet(uword &n, int &nbr_in, uword &var_in, vec &betaA, uvec &A, mat &x, mat &xt, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, double &lambda2, vec &xbar, sp_mat &spS, bool &usechol, uword &fun) ;

void add_var_enet(uword &n, int &nbr_in, uword &var_in, vec &betaA, uvec &A, sp_mat &x, sp_mat &xt, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, double &lambda2, vec &xbar, sp_mat &spS, bool &usechol, uword &fun) ;

void remove_var_enet(int &nbr_in, uvec &are_in, vec &betaA, uvec &A, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R,  uvec &null, bool &usechol, uword &fun) ;

template <typename any_mat>
void standardize(any_mat &x, vec &y, bool &intercept, bool &normalize, vec &penscale,
		 vec &xty, vec &normx, double &normy, vec &xbar, double &ybar) {

  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = trans(rowvec(mean(x, 0)));
    ybar = mean(y) ;
  } else {
    xbar = zeros(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(trans(sum(square(x),0)) - n * square(xbar));
    for (int i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
    xbar /= normx ;
  } else {
    normx = ones(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(penscale != 1)) {
    for (int i=0; i<n; i++) {
       x.row(i) /= penscale ;
    }
    xbar /= penscale;
  }

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
    for (int i=0;i<p;i++) {
       xty(i) -=  sum(y-ybar) * xbar(i);
    }
  } else {
    xty = trans(y.t()*x) ;
  }
}

#endif

