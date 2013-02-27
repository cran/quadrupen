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

#ifndef _quadrupen_UTILS_H
#define _quadrupen_UTILS_H
#include "utils.h"
#endif 

#include "quadratic.h"
#include "first_order.h"

RcppExport SEXP bounded_reg(SEXP X, SEXP XTY, SEXP STRUCT, SEXP LAMBDA1, SEXP LAMBDA2, SEXP XBAR, SEXP NORMX, SEXP WEIGHTS, SEXP NAIVE, SEXP EPS, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP BULLETPROOF) ;

#endif
