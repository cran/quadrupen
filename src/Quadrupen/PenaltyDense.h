#pragma once

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

enum class DenseNorm {L2, LINF};

using arma::vec;

template <DenseNorm norm> class DensePenalty {
public: 
  
  double gamma_ = 0 ; // Optional factor for non-concave penalties
  DensePenalty() {} ;

  vec    elt_norm  (const vec& x, const vec& w, double lambda=0) ;
  vec    elt_dual_norm  (const vec& x, const vec& w, double lambda=0) ;
  double pen_norm  (const vec& x, const vec& w, double lambda=0) ;
  double dual_norm (const vec& x, const vec& w, double lambda=0) ;
  vec proximal(const vec& x, double lambda, const vec& w) ;
  double lambda_max (const vec& XTy, const vec& w) ;
  vec optimality(const vec& x, double lambda, const vec& w)  ;
  
};

template<DenseNorm norm>
vec DensePenalty<norm>::optimality(const vec& grad, double lambda, const vec& w)  {
  return(elt_dual_norm(grad, w) - lambda) ;
}

template<DenseNorm norm>
double DensePenalty<norm>::lambda_max(const vec& XTy, const vec& w)  {
  return(dual_norm(XTy, w)) ;
}
