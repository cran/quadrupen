/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "Optimizer.h"
#include "PenaltyGroup.h"
#include "ActiveSetGroup.h"

using arma::vec;
using arma::mat;
using arma::uvec;
using arma::uword;
using arma::ones;

template <typename matrix, GroupSparseNorm norm>
class GroupOptimizer : public Optimizer {

public:

  GroupOptimizer() {} ;
  GroupOptimizer(GroupPenalty<norm>&, const List&) ;

  GroupPenalty<norm> penalty_ ;
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
      ActiveSetGroup<matrix>& set
  ) ;

  uword quadratic(
      vec &beta,
      const double lambda,
      const vec& weights,
      const vec &XTy,
      ActiveSetGroup<matrix> &set,
      const double& accuracy)  ;

};

template <typename matrix, GroupSparseNorm norm>
GroupOptimizer<matrix, norm>::GroupOptimizer(
    GroupPenalty<norm>& penalty, const List& control) :
  Optimizer(control) {
  penalty_  = penalty ;
}

template <typename matrix, GroupSparseNorm norm>
uword GroupOptimizer<matrix, norm>::quadratic(
    vec &beta,
    const double lambda,
    const vec& weights,
    const vec &XTy,
    ActiveSetGroup<matrix> &set,
    const double& tol) {

  double l1 = lambda * penalty_.alpha_ ;
  double l2 = lambda * (1-penalty_.alpha_) ;

  uword iter = 0;
  bool stable = false;
  const uword nb_active_groups = set.size_grp();

  while (!stable && iter < 15) {
    iter++;
    vec beta_old = beta;
    uword offset = 0; // go through active variables
    std::vector<uword> groups_to_remove_vec;

    for (uword k = 0; k < nb_active_groups; ++k) {
      uword sz   = set.grp_sizes_(set.G_(k)); // group size
      uvec ind_g = regspace<uvec>(offset, offset + sz - 1);

      vec beta_g = beta.subvec(offset, offset + sz - 1);

      // Get current group residuals:  XTy_g - X'X_g,A * beta_A
      vec res_g = XTy(ind_g) - set.XATXA_.rows(ind_g) * beta;
      res_g += set.XATXA_(ind_g, ind_g) * beta_g;

      if (norm == GroupSparseNorm::COOP) { // COOPERATIVE LASSO
        uvec idx_pos = find(beta_g >= 0);
        uvec idx_neg = find(beta_g < 0);

        mat Hg = set.XATXA_(ind_g, ind_g);
        if (!idx_pos.is_empty()) {
          double n_pos = arma::norm(beta_g(idx_pos), 2);
          for(uword i : idx_pos) Hg(i,i) += (l2 * weights(k)) / (n_pos + 1e-15);
        }
        if (!idx_neg.is_empty()) {
          double n_neg = arma::norm(beta_g(idx_neg), 2);
          for(uword i : idx_neg) Hg(i,i) += (l2 * weights(k)) / (n_neg + 1e-15);
        }

        arma::solve(beta_g, Hg, res_g, arma::solve_opts::fast);
      } else { // GROUP LASSO
        // Use cached Eigen decomposition
        double n_g = std::sqrt(arma::dot(beta_g, beta_g) + 1e-12);
        double mu = (l2 * weights(k)) / n_g;

        // Solve with
        vec D_mu_inv = 1.0 / (set.D_[k] + mu);
        beta_g = set.V_[k] * (D_mu_inv % (set.V_[k].t() * res_g));
      }

      // Soft-thresholding for Sparse Group
      if (l1 > 0) {
        for (uword i = 0; i < sz; ++i) {
          // Approximation of the curvature with diag of XTX
          double h_ii = set.XATXA_(ind_g(i), ind_g(i));
          beta_g(i) = arma::sign(beta_g(i)) * std::max(0.0, std::abs(beta_g(i)) - (l1 * weights(k)) / (h_ii + 1e-15));
        }
      }

      // Group deletion if applicable
      if (arma::norm(beta_g, 2) < 1e-10) {
        beta_g.zeros();
        groups_to_remove_vec.push_back(k);
      }

      beta.subvec(offset, offset + sz - 1) = beta_g;

      offset += sz; // go to next group
    }

    if (!groups_to_remove_vec.empty()) {
      uvec groups_to_remove = arma::conv_to<uvec>::from(groups_to_remove_vec);
      if (verbosity_) set.G_(groups_to_remove).print("\tremoving group") ;
      set.del_groups(groups_to_remove, beta);
      // Let the working set algorithm handle new active set (beta size now has changed...)
      return iter;
    }

    if (arma::norm(beta - beta_old, "inf") < tol) stable = true;

  }
  return iter;
}

template <typename matrix, GroupSparseNorm norm>
uword GroupOptimizer<matrix,norm>::working_set(
    vec& beta,
    vec& grad,
    const double& lambda,
    const vec& weights,
    const double& gamma,
    RegressionData<matrix> &data,
    ActiveSetGroup<matrix>& set) {

  if (verbosity_) Rprintf("\n current penalty = %f",lambda) ;
  if (verbosity_) Rprintf("\n nb active groups = %i\n", set.size_grp()) ;

  vec optimality = penalty_.optimality(grad, lambda, set.grp_sizes_, weights) ;
  uword grp_in = optimality.index_max() ; // highest violation of KKT conditions
  uword status = 0 ; iter_ = 0 ; bool success = true ;
  gap_ = std::max(0.0, optimality(grp_in)) ;
  J_ = arma::datum::inf ; D_ = arma::datum::inf ;

  double cached_L = -1.0 ; // Lipschitz constant cache; -1 means stale/not yet computed
  bool set_changed = true ; // active set changed since last Lipschitz computation

  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    double current_tol = 1e-7;

    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_grp_in_[grp_in] == 0 && optimality(grp_in) > 0) { // Is var_in already in the active set?
      set.add_group(grp_in, data) ;
      beta.insert_rows(beta.n_elem, set.grp_sizes_(grp_in)); // update the vector of active parameters
      if (algorithm_ ==  SolverType::QUADRA) {
        beta.tail(set.grp_sizes_(grp_in)).fill(- 1e-3 * arma::sign(grad(grp_in)));
      } else {
        beta.tail(set.grp_sizes_(grp_in)).fill(0.0);
      }
      if (verbosity_) {Rprintf("\tnewly added group %i\n",grp_in);}
      set_changed = true ;
    } else {
      set_changed = false ;
    }

    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ == SolverType::QUADRA) {
      inner_iter_.push_back(
        quadratic(beta, lambda, weights(set.G_), data.XTy_(set.A_), set, 1e-4)
      );
      grad = - data.XTy_ + set.XTXA_ * beta ;
    } else {
      if (set_changed) cached_L = estimate_lipschitz(set.XATXA_) ;
      auto prox = [this, &set, &weights](const vec& x, const double l) {
        return(penalty_.proximal(x, l, set.grp_sizes_(set.G_), weights(set.G_)));
      } ;
      vec beta_old = beta ;
      if (algorithm_ == SolverType::FISTA) {
        inner_iter_.push_back(
          fista(beta, lambda, data.XTy_(set.A_), set.XATXA_, prox, current_tol, 3000, cached_L)
        );
      } else if (algorithm_ == SolverType::PGD) {
        inner_iter_.push_back(
          pgd(beta, lambda, data.XTy_(set.A_), set.XATXA_, prox, current_tol, 3000, 3, cached_L)
        );
      }
      grad += set.XTXA_ * (beta - beta_old);
    }

    // VARIABLE DELETION IF APPLICABLE
    uvec vanish = find(
      penalty_.optimality(grad(set.A_), lambda, set.grp_sizes_(set.G_), weights(set.G_)) <= accuracy_ &&
        penalty_.elt_norm(beta, set.grp_sizes_(set.G_), ones(set.size_grp())) <= accuracy_/10
    ) ;
    if (!vanish.is_empty()) {
      if (verbosity_) set.G_(vanish).print("\tremoved group %i\n") ;
      set.del_groups(vanish, beta) ;
      set_changed = true ;
    }

    // OPTIMALITY TESTING
    optimality = penalty_.optimality(grad, lambda, set.grp_sizes_, weights) ;
    grp_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(grp_in)) ;

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
