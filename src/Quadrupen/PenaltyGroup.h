#pragma once

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>
#include "PenaltyUtils.h"

enum class GroupSparseNorm {L1L2, L1LINF, COOP};

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::zeros;

template <GroupSparseNorm norm> class GroupPenalty {
public: 
  
  double alpha_ = 0 ; // Optional sparse factor 
  GroupPenalty() {} ;
  GroupPenalty(double alpha) : alpha_(alpha){} ;
  
  double grp_norm  (const vec& x) ;
  double grp_norm_dual (const vec& x) ;
  vec    elt_norm  (const vec& x, const uvec& pk, const vec& wk) ;
  double pen_norm  (const vec& x, const uvec& pk, const vec& wk) ;
  double lambda_max (const vec& XTy, const uvec& pk, const vec& wk) ;
  vec proximal(const vec& x, double lambda, const uvec& pk, const vec& wk)  ;
  vec optimality(const vec& x, double lambda, const uvec& pk, const vec& wk)  ;
  
};

template<GroupSparseNorm norm>
vec GroupPenalty<norm>::elt_norm(const vec& x, const uvec& pk, const vec& wk) {
  
  vec  res = zeros<vec> (pk.n_elem) ;
  uword ind = 0 ;
  
  for (uword k=0; k<pk.n_elem; k++) {
    if (alpha_ > 0.0) res(k) += alpha_ * arma::norm(x.subvec(ind, ind + pk(k) - 1), 1);
    res(k) += (1 - alpha_) * wk(k) * grp_norm(x.subvec(ind, ind + pk(k) - 1));
    ind += pk(k);
  }
  
  return(res);
}

template<GroupSparseNorm norm>
double GroupPenalty<norm>::pen_norm(const vec& x, const uvec& pk, const vec& wk) {
  return(sum(elt_norm(x, pk, wk)));
}

template<GroupSparseNorm norm>
vec GroupPenalty<norm>::optimality(const vec& x, double lambda, const uvec& pk, const vec& wk)  {
  
  vec x_st = x;
  if (alpha_ > 0.0)
    x_st = soft_threshold(x_st, lambda * alpha_) ;
  
  vec res = zeros<vec>(pk.n_elem);
  uword ind = 0;
  for (uword k = 0; k < pk.n_elem; k++) {
    uword size = pk(k);
    res(k) = grp_norm_dual(x_st.subvec(ind, ind + size - 1));
    ind += size;
  }
  
  return(res / (wk  * (1 - alpha_)) - lambda) ;
  
}

template<GroupSparseNorm norm>
double GroupPenalty<norm>::lambda_max(const vec& XTy, const uvec& pk, const vec& wk) {
  double l_max = 0;
  
  // Look for maximal violation
  uword ind = 0 ; // index to go through the groups    
  for (uword k = 0; k < pk.n_elem; ++k) {
    uword size = pk(k);
    vec grad_k = XTy.subvec(ind, ind + size - 1);
    
    auto violation = [&](double l) {
      vec st = sign(grad_k) % max(abs(grad_k) - l * alpha_, zeros<vec>(grad_k.n_elem));
      return grp_norm_dual(st) - l * (1.0 - alpha_) * wk(k);
    };
    
    double low = 0, high = arma::norm(grad_k, "inf") / std::max(alpha_, 1e-3);
    high = std::max(high, grp_norm_dual(grad_k) / ((1.0-alpha_)*wk(k)));

    // Use bissection to find matching value
    for(int i=0; i<30; ++i) {
      double mid = (low + high) / 2.0;
      if (violation(mid) > 0) low = mid; 
      else high = mid;
    }

    if (high > l_max) l_max = high;
    
    ind += size;
  }
  
  return l_max ;
}
