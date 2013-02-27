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

RcppExport SEXP elastic_net(SEXP BETA0    ,
			    SEXP X        ,
			    SEXP Y        ,
			    SEXP STRUCT   ,
			    SEXP LAMBDA1  ,
			    SEXP N_LAMBDA ,
			    SEXP MIN_RATIO,
			    SEXP PENSCALE ,
			    SEXP LAMBDA2  ,
			    SEXP INTERCEPT,
			    SEXP NORMALIZE,
			    SEXP WEIGHTS  ,
			    SEXP NAIVE    ,
			    SEXP EPS      ,
			    SEXP MAXITER  ,
			    SEXP MAXFEAT  ,
			    SEXP FUN      ,
			    SEXP VERBOSE  ,
			    SEXP SPARSE   ,
			    SEXP USECHOL  ,
			    SEXP MONITOR  ) ;

#endif
