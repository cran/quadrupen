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

  // Exact minimizer over one active group of 1/2 b'Hb - r'b + penalty(b), with H = V diag(d) V'
  vec block_solve(
      const vec& r,
      const mat& H,
      const vec& d,
      const mat& V,
      const vec& beta_g,
      const double lambda,
      const double w,
      const double& accuracy) ;

};

template <typename matrix, GroupSparseNorm norm>
GroupOptimizer<matrix, norm>::GroupOptimizer(
    GroupPenalty<norm>& penalty, const List& control) :
  Optimizer(control) {
  penalty_  = penalty ;
}

template <typename matrix, GroupSparseNorm norm>
vec GroupOptimizer<matrix, norm>::block_solve(
    const vec& r,
    const mat& H,
    const vec& d,
    const mat& V,
    const vec& beta_g,
    const double lambda,
    const double w,
    const double& accuracy) {

  uword sz = r.n_elem ;
  uvec pk = {sz} ;
  vec wk = {w} ;

  // b = 0 is optimal iff the (soft-thresholded) partial residual lies in the dual ball
  if (penalty_.optimality(r, lambda, pk, wk)(0) <= 0.0) return zeros<vec>(sz) ;

  const double mu = lambda * (1.0 - penalty_.alpha_) * w ;
  const vec dd = arma::clamp(d, 0.0, arma::datum::inf) ;

  if (norm == GroupSparseNorm::L1L2 && penalty_.alpha_ == 0.0 && dd.max() > 0.0) {
    // Group-Lasso: b = (H + mu/t I)^{-1} r with t = ||b||, i.e. t solves the secular equation
    //   g(t) = sum_i c_i^2 / (d_i t + mu)^2 - 1 = 0,  c = V'r.
    // g is convex and decreasing with g(0) > 0 (b != 0), so Newton iterations started on the
    // left of the root increase monotonically towards it.
    const vec c2 = arma::square(V.t() * r) ;
    auto g = [&](double t) { return accu(c2 / arma::square(dd * t + mu)) - 1.0 ; } ;
    double t = arma::norm(beta_g, 2) ;
    if (!(t > 0.0) || g(t) < 0.0) t = 0.0 ;
    for (uword it = 0; it < 100; ++it) {
      const vec den = dd * t + mu ;
      const double gt  = accu(c2 / arma::square(den)) - 1.0 ;
      const double dgt = -2.0 * accu(c2 % dd / arma::pow(den, 3)) ;
      if (dgt >= 0.0) break ;
      const double step = - gt / dgt ;
      t += step ;
      if (std::abs(step) <= 1e-12 * t) break ;
    }
    return V * ((V.t() * r) % (t / (dd * t + mu))) ;
  }

  // Other penalties (sparse-group, coop, l1/linf): FISTA with restart on the block,
  // with the exact proximal operator of the group penalty
  const double L = std::max(dd.max(), 1e-12) ;
  vec b = beta_g, y = beta_g, b_old ;
  double t0 = 1.0 ;
  for (uword it = 0; it < 10000; ++it) {
    b_old = b ;
    b = penalty_.proximal(y - (H * y - r) / L, lambda / L, pk, wk) ;
    if (L * arma::norm(b - b_old, 2) <= accuracy) break ;
    if (dot(y - b, b - b_old) > 0) {
      t0 = 1.0 ; y = b ;
    } else {
      double t1 = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t0 * t0)) ;
      y = b + ((t0 - 1.0) / t1) * (b - b_old) ;
      t0 = t1 ;
    }
  }
  return b ;
}

template <typename matrix, GroupSparseNorm norm>
uword GroupOptimizer<matrix, norm>::quadratic(
    vec &beta,
    const double lambda,
    const vec& weights,
    const vec &XTy,
    ActiveSetGroup<matrix> &set,
    const double& tol) {

  // Block coordinate descent over the active groups, with an exact solve for each block
  const uword nb_active_groups = set.size_grp();
  if (nb_active_groups == 0) return 0 ;

  uvec first(nb_active_groups), last(nb_active_groups) ;
  uword offset = 0 ;
  for (uword k = 0; k < nb_active_groups; ++k) {
    first(k) = offset ;
    offset  += set.grp_sizes_(set.G_(k)) ;
    last(k)  = offset - 1 ;
  }

  // For the Group-Lasso, the penalty is smooth away from zero: each sweep is followed by a
  // damped Newton step on the nonzero groups, which fixes the slow (linear) convergence of
  // block coordinate descent when groups are strongly correlated
  const bool newton = (norm == GroupSparseNorm::L1L2 && penalty_.alpha_ == 0.0) ;

  uword sweep = 0 ;
  double max_delta = arma::datum::inf ;
  while (max_delta > tol * std::max(1.0, arma::norm(beta, "inf")) && sweep < 1000) {
    sweep++ ;
    max_delta = 0.0 ;
    for (uword k = 0; k < nb_active_groups; ++k) {
      const vec beta_g = beta.subvec(first(k), last(k)) ;
      const mat H = set.XATXA_.submat(first(k), first(k), last(k), last(k)) ;

      // partial residual of the group: X'y_g - X'X_{g,A} beta + H beta_g
      const vec r = XTy.subvec(first(k), last(k)) - set.XATXA_.rows(first(k), last(k)) * beta + H * beta_g ;

      vec b ;
      if (set.use_evd_) {
        b = block_solve(r, H, set.D_[k], set.V_[k], beta_g, lambda, weights(k), tol) ;
      } else {
        vec d ; mat V ;
        eig_sym(d, V, H) ;
        b = block_solve(r, H, d, V, beta_g, lambda, weights(k), tol) ;
      }
      max_delta = std::max(max_delta, arma::norm(b - beta_g, "inf")) ;
      beta.subvec(first(k), last(k)) = b ;
    }

    if (newton) {
      // Variables of the nonzero groups, and group boundaries within them
      std::vector<uword> idx_vec, gfirst, glast ;
      std::vector<double> gmu ;
      for (uword k = 0; k < nb_active_groups; ++k) {
        if (!arma::any(beta.subvec(first(k), last(k)))) continue ;
        gfirst.push_back(idx_vec.size()) ;
        for (uword i = first(k); i <= last(k); ++i) idx_vec.push_back(i) ;
        glast.push_back(idx_vec.size() - 1) ;
        gmu.push_back(lambda * weights(k)) ;
      }
      if (!idx_vec.empty()) {
        const uvec idx = arma::conv_to<uvec>::from(idx_vec) ;
        const mat Hn = set.XATXA_(idx, idx) ;
        const vec rn = XTy(idx) ;
        vec b = beta(idx) ;

        auto objective = [&](const vec& x) {
          double val = 0.5 * dot(x, Hn * x) - dot(rn, x) ;
          for (uword g = 0; g < gfirst.size(); ++g) val += gmu[g] * arma::norm(x.subvec(gfirst[g], glast[g]), 2) ;
          return val ;
        } ;

        vec grad_n = Hn * b - rn ;
        mat K = Hn ;
        for (uword g = 0; g < gfirst.size(); ++g) {
          const vec bg = b.subvec(gfirst[g], glast[g]) ;
          const double t = arma::norm(bg, 2) ;
          const vec u = bg / t ;
          grad_n.subvec(gfirst[g], glast[g]) += gmu[g] * u ;
          K.submat(gfirst[g], gfirst[g], glast[g], glast[g]) +=
            (gmu[g] / t) * (arma::eye(bg.n_elem, bg.n_elem) - u * u.t()) ;
        }
        vec dir ;
        if (arma::solve(dir, K, -grad_n, arma::solve_opts::likely_sympd + arma::solve_opts::fast)) {
          const double slope = dot(grad_n, dir) ;
          if (slope < 0.0) {
            // Armijo backtracking line search
            const double f0 = objective(b) ;
            double step = 1.0 ;
            while (objective(b + step * dir) > f0 + 1e-4 * step * slope && step > 1e-10) step *= 0.5 ;
            if (step > 1e-10) {
              beta(idx) = b + step * dir ;
              max_delta = std::max(max_delta, step * arma::norm(dir, "inf")) ;
            }
          }
        }
      }
    }
    if (sweep % 10 == 0) R_CheckUserInterrupt() ;
  }

  // Remove the groups set exactly to zero
  std::vector<uword> zero_groups ;
  for (uword k = 0; k < nb_active_groups; ++k) {
    if (!arma::any(beta.subvec(first(k), last(k)))) zero_groups.push_back(k) ;
  }
  if (!zero_groups.empty()) {
    uvec groups_to_remove = arma::conv_to<uvec>::from(zero_groups) ;
    if (verbosity_) set.G_(groups_to_remove).t().print("\tremoving groups") ;
    set.del_groups(groups_to_remove, beta) ;
  }

  return sweep ;
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
  uword status = 0 ; iter_ = 0 ; bool success = true ;
  gap_ = std::max(0.0, optimality.max()) ;
  J_ = arma::datum::inf ; D_ = arma::datum::inf ;

  double cached_L = -1.0 ; // Lipschitz constant cache; -1 means stale/not yet computed
  bool set_changed = true ; // active set changed since last Lipschitz computation

  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    double current_tol = 1e-7;

    // GROUP ACTIVATION IF APPLICABLE: the largest KKT violators among inactive groups,
    // stopping once more than maxfeat variables are active (which stops the path)
    uvec grps_in = select_violators(optimality, set.is_grp_in_, accuracy_, max_add_) ;
    for (uword grp_in : grps_in) {
      if (set.size() > maxfeat_) break ;
      set.add_group(grp_in, data) ;
      beta = arma::join_cols(beta, arma::zeros<vec>(set.grp_sizes_(grp_in))) ;
      if (verbosity_) {Rprintf("\tnewly added group %i\n",grp_in);}
      set_changed = true ;
    }

    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ == SolverType::QUADRA) {
      inner_iter_.push_back(
        quadratic(beta, lambda, weights(set.G_), data.XTy_(set.A_), set, 1e-9)
      );
      grad = - data.XTy_ + set.XTXA_times(beta) ;
    } else {
      if (set_changed) { cached_L = estimate_lipschitz(set.XATXA_) ; set_changed = false ; }
      const uvec pkA = set.grp_sizes_(set.G_) ; // computed once, not at each inner iteration
      const vec  wkA = weights(set.G_) ;
      const vec XTyA = data.XTy_(set.A_) ;
      auto prox = [this, &pkA, &wkA](const vec& x, const double l) {
        return(penalty_.proximal(x, l, pkA, wkA));
      } ;
      vec beta_old = beta ;
      if (algorithm_ == SolverType::FISTA) {
        inner_iter_.push_back(
          fista(beta, lambda, XTyA, set.XATXA_, prox, current_tol, 3000, cached_L)
        );
      } else if (algorithm_ == SolverType::PGD) {
        inner_iter_.push_back(
          pgd(beta, lambda, XTyA, set.XATXA_, prox, current_tol, 3000, 3, cached_L)
        );
      }
      grad += set.XTXA_times(beta - beta_old);
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
    gap_ = std::max(0.0, optimality.max()) ;

    if (monitoring_ > 0) {
      optimality_violation(beta, grad, lambda, gamma, data.XTy_(set.A_), set.XATXA_, data.norm_y_, set.A_, monitoring_) ;
      J_vec_.push_back(J_) ;
      D_vec_.push_back(D_) ;
    }

  }
  if (verbosity_) Rprintf("\tcurrent gap = %f\n",gap_) ;

  // Checking convergence status
  if (gap_ > accuracy_)      { status = 1 ; }
  if (set.size() > maxfeat_) { status = 2 ; }
  if (!success)              { status = 3 ; }

  return status ;
}
