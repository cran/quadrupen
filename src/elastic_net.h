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
#include "coordinate.h"
#include "proximal.h"

RcppExport SEXP elastic_net(SEXP X, SEXP XTY, SEXP S, SEXP LAMBDA, SEXP PENSCALE, SEXP XBAR, SEXP YBAR, SEXP NORMX, SEXP NORMY, SEXP WEIGHTS, SEXP NAIVE, SEXP EPS, SEXP ZERO, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP SPARSE, SEXP MONITOR) ;

#endif
