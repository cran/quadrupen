/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "Optimizer.h"
#include "PenaltySparse.h"
#include "ActiveSet.h"

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::mat;
using arma::regspace;

template <typename matrix, SparseNorm norm>
class SparseOptimizer: public Optimizer {

public:

  SparseOptimizer() {} ;
  SparseOptimizer(SparsePenalty<norm>&, const List&) ;

  SparsePenalty<norm> penalty_ ;
  using Optimizer::algorithm_  ;
  using Optimizer::accuracy_   ;
  using Optimizer::maxiter_    ;
  using Optimizer::maxfeat_    ;
  using Optimizer::verbosity_  ;
  using Optimizer::monitoring_ ;
  using Optimizer::iter_       ;
  using Optimizer::inner_iter_ ;
  using Optimizer::gap_        ;
  using Optimizer::J_          ;
  using Optimizer::D_          ;
  using Optimizer::J_vec_      ;
  using Optimizer::D_vec_      ;
  using Optimizer::optimality_violation ;
  using Optimizer::fista ;
  using Optimizer::pgd ;

  uword working_set(
      vec& beta,
      vec& grad,
      const double& lambda,
      const vec& weights,
      const double& gamma,
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set
  ) ;

  uword quadratic(
      vec &beta,
      const double &lambda,
      const vec &weights,
      const vec &XTy,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) ;

};

template <typename matrix, SparseNorm norm>
SparseOptimizer<matrix, norm>::SparseOptimizer(
    SparsePenalty<norm>& penalty, const List& control) :
    Optimizer(control) {
    penalty_ = penalty ;
  }

template <typename matrix, SparseNorm norm>
uword SparseOptimizer<matrix, norm>::quadratic(
    vec &beta,
    const double &lambda,
    const vec &weights,
    const vec &XTy,
    ActiveSet<matrix> &set,
    const double& accuracy,
    const uword& max_iter) {

  uword iter = 0 ;
  bool convergence = false; // check for sign stability and LLA consistency

  while (!convergence && iter < 50 && set.size() > 0) { // Max 10 swaps
    iter++;
    // Local weights for local linear approximation (LLA)
    // For MCP/SCAD, effective weights changes according to the current beta
    vec local_w = penalty_.derivative(beta, lambda, weights.elem(set.A_));
    vec theta = arma::sign(beta);
    theta.replace(0.0, 1.0);  // Handle variable that has just been zeroed

    // Solving the quadratic problem (KKT Newton)
    // (XA'XA) beta = XA'y - local_w * sign(beta)
    vec beta_new; // candidate for next step
    if (set.use_chol_) {
      vec rhs = XTy(set.A_) - local_w % theta;
      // Step 1 - Forward Substitution
      vec tmp = arma::solve(trimatl(set.R_.t()), rhs);
      // Step 2 - Backward Substitution
      beta_new = arma::solve(trimatu(set.R_), tmp);
    } else {
      beta_new = beta ; // warm start for CG
      this->conjugate_gradient(beta_new, set.XATXA_,
                               XTy(set.A_) - local_w % theta,
                               accuracy, max_iter);
    }

    // Check for swapping variables / sign stability
    uvec swap_idx = find(arma::sign(beta_new) != theta && abs(beta) > accuracy);
    if (swap_idx.is_empty()) { // No swap: check for convergence
      double diff = arma::norm(beta_new - beta, 2); // for MCP and SCAD
      beta = beta_new;
      if (diff < accuracy) convergence = true;
    } else {
      // Find the first variable hitting zero and its interpolating ratio
      vec ratios = -beta(swap_idx) / (beta_new(swap_idx) - beta(swap_idx));
      uword i_min = ratios.index_min();
      uword idx_to_remove = swap_idx[i_min];

      // Interpolate by moving all variables to this point
      beta = beta + ratios(i_min) * (beta_new - beta);
      beta[idx_to_remove] = 0.0;

      // Remove the incriminated variable
      if (verbosity_) Rprintf("\tremoving variables %i\n", set.A_(idx_to_remove)) ;
      set.del_var(idx_to_remove, beta) ;
    }
  }

  return(iter) ;

}

template <typename matrix, SparseNorm norm>
uword SparseOptimizer<matrix,norm>::working_set(
    vec& beta,
    vec& grad,
    const double& lambda,
    const vec& weights,
    const double& gamma,
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set) {

  if (verbosity_) Rprintf("\n current penalty = %f",lambda) ;
  if (verbosity_) Rprintf("\n nb active variables = %i\n", set.size()) ;

  vec optimality = penalty_.optimality(grad, lambda, weights, beta, set.A_);
  uword var_in = optimality.index_max() ; // highest violation of KKT conditions
  uword status = 0 ; iter_ = 0 ; bool success = true ;
  gap_ = std::max(0.0, optimality(var_in)) ;
  J_ = arma::datum::inf ; D_ = arma::datum::inf ;

  double cached_L = -1.0 ; // Lipschitz constant cache; -1 means stale/not yet computed
  bool set_changed = true ; // active set changed since last Lipschitz computation

  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;

    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_in_[var_in] == 0) { // Is var_in already in the active set?
      set.add_var(var_in, data) ;
      beta.insert_rows(beta.n_elem, 1); // update the vector of active parameters
      if (algorithm_ ==  SolverType::QUADRA) {
        beta.tail(1).fill(- 1e-3 * arma::sign(grad(var_in)));
      } else {
        beta.tail(1).fill(0.0);
      }
      if (verbosity_) {Rprintf("\tnewly added variable %i\n",var_in);}
      set_changed = true ;
    } else {
      set_changed = false ;
    }

    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ ==  SolverType::QUADRA) { // Newton-based solver
      inner_iter_.push_back(
        quadratic(beta, lambda, weights, data.XTy_, set, 1e-9, 1000)
      );
      grad = - data.XTy_ + set.XTXA_ * beta ;
    }
    else { // Proximal-based solvers
      if (set_changed) cached_L = estimate_lipschitz(set.XATXA_) ;
      auto prox = [this, &set, &weights](const vec& x, double l) {
        return(penalty_.proximal(x, l, weights.elem(set.A_)));
      };
      vec beta_old = beta ;
      if (algorithm_ == SolverType::FISTA) {
        inner_iter_.push_back(
          fista(beta, lambda, data.XTy_.elem(set.A_), set.XATXA_, prox, 1e-7, 3000, cached_L)
        );
      } else if (algorithm_ == SolverType::PGD) {
        inner_iter_.push_back(
          pgd(beta, lambda, data.XTy_.elem(set.A_), set.XATXA_, prox, 1e-7, 3000, 5, cached_L)
        );
      }
      grad += set.XTXA_ * (beta - beta_old); // Incremental update of the gradient

      uvec local_A = regspace<uvec>(0, set.size() - 1); // local indices
      vec kkt_res  = penalty_.optimality(
        grad.elem(set.A_), lambda, weights.elem(set.A_),
        beta, local_A
      );
      uvec vanish = find(kkt_res <= accuracy_ &&
          abs(beta) < accuracy_/10 * weights.elem(set.A_)
        );
      if (!vanish.is_empty()) {
        if (verbosity_) {set.A_(vanish).t().print("Removing variables");}
        set.del_vars(vanish, beta) ;
        set_changed = true ;
      }
    }

    // OPTIMALITY TESTING
    optimality = penalty_.optimality(grad, lambda, weights, beta, set.A_);
    var_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(var_in)) ;

    if (monitoring_ > 0) {
      optimality_violation(beta, grad, lambda, gamma, data.XTy_(set.A_), set.XATXA_, data.norm_y_, set.A_, monitoring_) ;
      J_vec_.push_back(J_) ;
      D_vec_.push_back(D_) ;
    }

  }
  if (verbosity_) Rprintf("\tcurrent gap = %f\n",gap_) ;

  // Checking convergence status
  if (iter_ >= maxiter_)     { status = 1 ; }
  if (set.size() > maxfeat_) { status = 2 ; }
  if (!success)              { status = 3 ; }

  return status ;
}
