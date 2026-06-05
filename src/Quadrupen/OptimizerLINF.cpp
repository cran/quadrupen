/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "OptimizerLINF.h"

using namespace Rcpp;
using namespace arma;

OptimizerLINF::OptimizerLINF(
  DensePenalty<DenseNorm::LINF>& penalty, const List& control) : 
  Optimizer(control) 
{penalty_ = penalty ;}

uword OptimizerLINF::quadratic_breg(
  vec& beta,
  vec &grad,
  const double& lambda,
  const vec& weights,
  RegressionData<mat> &data,
  uvec& U, // unbounded variables
  const double& accuracy,
  const uword& max_iter) {
  
  uvec B = regspace<uvec>(0,beta.n_elem-1) ; B.shed_rows(U) ;
  grad = -data.XTy_ + data.XTX_ * beta ;
  vec theta = -sign(grad(B)) ; // sign of the guys on the boundary
  
  uword iter = 0, iter_in= 0 ; // count the number of systems solved
  bool success = false ;
  const uword max_boundary_iter = beta.n_elem + 1; // at most p variables can cross
  while ((!success) && (iter < max_boundary_iter)) {
    
    iter++;
    
    // SOLVE THE QUADRATIC PROBLEM
    vec XX_B = data.XTX_.cols(B) * theta;
    if (U.is_empty()) {
      double b  = dot(theta, data.XTy_(B)) - lambda * mean(weights(B));
      beta(B) = theta * (b/sum(theta % XX_B(B),0)) ;
    } else {
      // Constructing the system (KKT)
      vec tmp = join_cols(XX_B(U),sum(theta % XX_B(B),0)) ;
      mat XX = join_rows(join_cols(data.XTX_(U, U), trans(XX_B(U))), tmp);
      vec b = data.XTy_(U) ; b.resize(b.n_elem + 1) ;
      b.tail(1) = dot(theta, data.XTy_(B)) - lambda * mean(weights(B));
      vec pen_arg = beta(B);
      double bound = max(abs(pen_arg)) ;
      tmp = join_cols(beta(U), ones(1) * bound);
      iter_in = this->conjugate_gradient(tmp, XX, b, accuracy, max_iter) ;
      beta(B) = theta * tmp.tail(1) ;
      beta(U) = tmp.subvec(0,tmp.n_elem-2) ;
    }
    
    // Handling guys reaching the boundary
    vec pen_arg = beta(B);
    double bound = max(abs(pen_arg)) ; // current boundary
    uvec ind_toB = find(abs(beta(U)) > bound) ;
    if (!ind_toB.is_empty()) {
      uvec toB = U(ind_toB) ;
      U.shed_rows(ind_toB) ;
      B = join_cols(B, toB);
      beta(B) = bound * sign(beta(B));
      theta = sign(beta(B));
    } else {
      success = true ;
    }
  }
  
  if (!success)
    throw std::runtime_error("Failed to resolve boundary transitions (too many simultaneous crossings).");

  // Guys leaving the boundary after optimization
  grad = -data.XTy_ + data.XTX_ * beta ;
  uvec ind_toA  = find(theta == sign(grad(B)));
  if (!ind_toA.is_empty()) {
    uvec toA = B(ind_toA) ;
    U.resize(U.n_elem + toA.n_elem);
    U.tail(toA.n_elem) = toA;
    B.shed_rows(ind_toA) ;
    if (B.is_empty()) {
      throw std::runtime_error("Every variable left the boudary. Try more regularization.");
    }
  }
  
  return(iter_in) ;
  
}
