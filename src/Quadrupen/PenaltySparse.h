#pragma once

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

enum class SparseNorm {L1, MCP, SCAD};

using arma::vec;
using arma::uvec;

template <SparseNorm norm> class SparsePenalty {
  public: 
    
    double eta_ = 0 ; // Optional factor for non-concave penalties
    SparsePenalty() {} ;
    SparsePenalty(double eta) : eta_(eta){} ;
    
    vec    elt_norm  (const vec& x, const vec& w, double lambda=0) ;
    vec    elt_dual_norm  (const vec& x, const vec& w, double lambda=0) ;
    double pen_norm  (const vec& x, const vec& w, double lambda=0) ;
    double dual_norm (const vec& x, const vec& w, double lambda=0) ;
    vec proximal(const vec& x, double lambda, const vec& w) ;
    double lambda_max (const vec& XTy, const vec& w) ;
    // vec optimality(const vec& x, double lambda, const vec& w)  ;
    vec optimality(const vec& grad, double lambda, const vec& w,
                   const vec& beta = vec(), const uvec& A = uvec()) const;
    
    vec derivative(const vec& beta, double lambda, const vec& w) ;

};

template<SparseNorm norm>
double SparsePenalty<norm>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(accu(elt_norm(x, w)));
}

template<SparseNorm norm>
vec SparsePenalty<norm>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) / w);
}

template<SparseNorm norm>
double SparsePenalty<norm>::dual_norm(const vec& x, const vec& w, double lambda) {
  return arma::max(elt_dual_norm(x, w));
}

template<SparseNorm norm>
  double SparsePenalty<norm>::lambda_max(const vec& XTy, const vec& w)  {
    return(dual_norm(XTy, w)) ;
  }

