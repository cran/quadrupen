/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_ELASTICNET_H
#define _quadrupen_ELASTICNET_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <sys/time.h>
#include <RcppArmadillo.h>

#include "utils.h"
#include "quadratic.h"
#include "first_order.h"

using namespace Rcpp;
using namespace arma;

void cholupdate(mat &R, mat& XtX) ;

void choldowndate(mat &R, int j) ;

void bound_to_optimal(vec &betaA, mat &xAtxA, vec &xty, vec &grd, double &lambda1, double &lambda2,
		      double &normy, uvec &A, int &monitor, vec &J_hat, vec &D_hat) ;

void add_var_enet(uword &n, int &nbr_in, uword &var_in, vec &betaA, uvec &A, mat &x, mat &xtxA, mat &xAtxA, mat &xtxw,
		  mat &R, double &lambda2, vec &xbar, vec &Xx, uvec &Xi, uvec &Xp, uvec &Xnp, uvec &j_nz,
		  vec &Sx, uvec &Si, uvec& Sp, bool &sparse, uword &fun) ;

void remove_var_enet(int &nbr_in, uvec &are_in, vec &betaA, uvec &A, mat &xtxA, mat &xAtxA, mat &xtxw,
		     mat &R,  uvec &null, uword &fun) ;

RcppExport SEXP elastic_net(SEXP X, SEXP XTY, SEXP S, SEXP LAMBDA1, SEXP LAMBDA2, SEXP XBAR, SEXP NORMX, SEXP NORMY, SEXP WEIGHTS, SEXP NAIVE, SEXP EPS, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP SPARSE, SEXP MONITOR) ;

#endif
