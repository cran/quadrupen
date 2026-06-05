/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "PenaltyGroup.h"
#include "PenaltyUtils.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// L1/L2 NORM A.K.A GROUP-LASSO
template<>
double GroupPenalty<GroupSparseNorm::L1L2>::grp_norm(const vec& x) {
  return(arma::norm(x, 2)) ;
}

template<>
double GroupPenalty<GroupSparseNorm::L1L2>::grp_norm_dual(const vec& x) {
  return(arma::norm(x, 2)) ;
}

template<>
vec GroupPenalty<GroupSparseNorm::L1L2>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {

  vec res = x;
  
  double l1_threshold = lambda * alpha_;
  double l2_grp_const = lambda * (1.0 - alpha_);
  
  if (l1_threshold > 0) {
    res = soft_threshold(x, l1_threshold);
  }
  
  uword ind = 0;
  for (uword k = 0; k < pk.n_elem; k++) {
    uword size = pk(k);
    auto group_slice = res.subvec(ind, ind + size - 1);
    
    double group_norm = grp_norm(group_slice);
    double threshold_g = l2_grp_const * wk(k);
    if (group_norm > threshold_g && group_norm > 0) {
      double shrink = 1.0 - (threshold_g / group_norm);
      group_slice *= shrink;
    } else {
      group_slice.zeros();
    }
    ind += size;
  }
  
  return res;
}

// ______________________________________________________
// L1/LINF NORM A.K.A GROUP-LASSO type 2
// 

template<>
double GroupPenalty<GroupSparseNorm::L1LINF>::grp_norm(const vec& x) {
  return(norm(x, "inf")) ;
}

template<>
double GroupPenalty<GroupSparseNorm::L1LINF>::grp_norm_dual(const vec& x) {
  return(norm(x, 1)) ;
}

template<>
vec GroupPenalty<GroupSparseNorm::L1LINF>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {

  vec res = x;
  
  double l1_threshold = lambda * alpha_;
  double l2_grp_const = lambda * (1.0 - alpha_);
  
  if (l1_threshold > 0) {
    res = soft_threshold(x, l1_threshold);
  }

  uword ind = 0;
  for (uword k = 0; k < pk.n_elem; k++) {
    uword size = pk(k);
    vec x_k = res.subvec(ind, ind + size - 1);
    double l2_w = l2_grp_const * wk(k);

    if (accu(abs(x_k)) > l2_grp_const * wk(k)) {
      // L1 ball projection (Moreau : prox_inf(v) = v - proj_L1(v))
      vec u = sort(arma::abs(x_k), "descend");
      vec css = cumsum(u);
      
      // Projection threshold rho
      double rho = 0;
      for (uword j = 0; j < size; ++j) {
        double val = (css(j) - l2_w) / (j + 1.0);
        if (j < size - 1 && u(j+1) <= val) {
          rho = val;
          break;
        }
        if (j == size - 1) rho = val;
      }
      
      // Proximal L_inf : x_k = sign(x_k) * min(|x_k|, rho)
      x_k = sign(x_k) % min(abs(x_k), rho * ones(size) );
    }
    
    res.subvec(ind, ind + size - 1) = x_k;
    ind += size;
  }
  return res;
}

// ______________________________________________________
// COOP(ERATIVE) NORM A.K.A COOPERATIVE-LASSO
// 

template<>
double GroupPenalty<GroupSparseNorm::COOP>::grp_norm(const vec& x) {
  double pos_sq = 0.0;
  double neg_sq = 0.0;
  
  for(double val : x) {
    if (val > 0) pos_sq += val * val;
    else neg_sq += val * val;
  }
  
  return std::sqrt(pos_sq) + std::sqrt(neg_sq);
}

template<>
double GroupPenalty<GroupSparseNorm::COOP>::grp_norm_dual(const vec& x) {
  double pos_sq = 0.0;
  double neg_sq = 0.0;
  
  for(double val : x) {
    if (val > 0) pos_sq += val * val;
    else neg_sq += val * val;
  }
  
  return( (pos_sq > neg_sq) ? std::sqrt(pos_sq) : std::sqrt(neg_sq) ) ;
}

template<>
vec GroupPenalty<GroupSparseNorm::COOP>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {
  
  vec res = x;
  
  double l1_threshold = lambda * alpha_;
  double l2_grp_const = lambda * (1.0 - alpha_);
  
  if (l1_threshold > 0) {
    res = soft_threshold(x, l1_threshold);
  }
  
  uword ind = 0;
  for (uword k=0; k<pk.n_elem; k++) {
    uword size = pk(k) ;
    vec x_k = res.subvec(ind, ind + size - 1);
    double l2_w = l2_grp_const * wk(k);
    
    double pos_sq = 0.0;
    double neg_sq = 0.0;

    for(double val : x_k) {
      if (val > 0) pos_sq += val * val;
      else neg_sq += val * val;
    }

    double norm_pos = std::sqrt(pos_sq);
    double norm_neg = std::sqrt(neg_sq);
    double shrink_pos = (norm_pos > l2_w) ? (1.0 - l2_w / norm_pos) : 0.0;
    double shrink_neg = (norm_neg > l2_w) ? (1.0 - l2_w / norm_neg) : 0.0;

    for (uword j=ind; j<(ind+size); j++) {
      if (res[j] > 0) res[j] *= shrink_pos ;
       else if (res[j] < 0) res[j] *= shrink_neg ;    
    }
    
    ind += size;
  }

  return(res);
  
}
