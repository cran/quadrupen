/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

// ====================================================
// Group-Sparse Regularizers

#pragma once

#include "SparsifyingRegularizer.h"
#include "OptimizerGroup.h"

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::sp_mat;
using arma::umat;
using arma::ones;
using Rcpp::List;
using Rcpp::as;
using std::vector;

template <typename matrix, GroupSparseNorm norm>
class GroupSparseRegularizer :
  public SparsifyingRegularizer<matrix> {
  public:

    using SparsifyingRegularizer<matrix>::intercept_ ;
    using SparsifyingRegularizer<matrix>::gamma_     ;
    using SparsifyingRegularizer<matrix>::data_      ;
    using SparsifyingRegularizer<matrix>::beta_      ;
    using SparsifyingRegularizer<matrix>::grad_      ;
    using SparsifyingRegularizer<matrix>::lambda_factor_ ;
    using SparsifyingRegularizer<matrix>::get_lambda_seq ;
    using SparsifyingRegularizer<matrix>::nzeros_    ;
    using SparsifyingRegularizer<matrix>::debiased_  ;
    using SparsifyingRegularizer<matrix>::active_    ;
    using SparsifyingRegularizer<matrix>::intercept_debiased_ ;
    using SparsifyingRegularizer<matrix>::beta_debiased_ ;
    using typename SparsifyingRegularizer<matrix>::StepResult ;
    using typename SparsifyingRegularizer<matrix>::Diagnostics ;

    ActiveSetGroup<matrix>            set_     ;
    GroupPenalty<norm>                penalty_ ;
    GroupOptimizer<matrix,norm>       solver_  ;

    double get_lambda_max() {
      return penalty_.lambda_max(data_.XTy_, set_.grp_sizes_, lambda_factor_) ;
    }

    GroupSparseRegularizer(RegressionData<matrix>, const uvec&, const List&, const List&) ;

    double get_df() override ;

    // ── SparsifyingRegularizer hooks ─────────────────────────────────────────

    ActiveSet<matrix>& current_set() override { return set_ ; }

    StepResult run_solver(double lambda) override {
      uword status = solver_.working_set(beta_, grad_, lambda, lambda_factor_, gamma_, data_, set_) ;
      return { status, solver_.gap_, solver_.iter_ } ;
    }

    Diagnostics solver_diagnostics() const override {
      return { solver_.inner_iter_, solver_.J_vec_, solver_.D_vec_ } ;
    }

} ;

template <typename matrix, GroupSparseNorm norm>
GroupSparseRegularizer<matrix,norm>::GroupSparseRegularizer(
  RegressionData<matrix> data, const uvec& group_ind, const List& regParam, const List& control) :
  SparsifyingRegularizer<matrix>::SparsifyingRegularizer(std::move(data), regParam) {

    // Scale the structuring matrix according to the amount of l2 penalty
    data_.scale_struct(gamma_) ;

    // Initialize the active set and the gradient
    grad_ = - data_.XTy_ ;
    set_ = ActiveSetGroup(data_, group_ind, as<bool>(control["factmat"])) ;

    // Set up the penalty
    penalty_ = GroupPenalty<norm>(as<double>(regParam["alpha"])) ;
    get_lambda_seq(penalty_.lambda_max(data_.XTy_, set_.grp_sizes_, lambda_factor_), regParam) ;

    // Set up the optimizer
    solver_ = GroupOptimizer<matrix,norm>(penalty_, control) ;
  }

template <typename matrix, GroupSparseNorm norm>
double GroupSparseRegularizer<matrix,norm>::get_df() {

  double df = data_.centered_ ;

  if (set_.size_grp() > 0) {
    // approximate degrees of freedom
    vec active_grp_norm     = penalty_.elt_norm(beta_,          set_.grp_sizes_(set_.G_), ones(set_.size_grp())) ;
    vec active_grp_norm_ols = penalty_.elt_norm(beta_debiased_, set_.grp_sizes_(set_.G_), ones(set_.size_grp())) / (1 + gamma_) ;

    df = df +
      accu(1 + (active_grp_norm / active_grp_norm_ols) % (set_.grp_sizes_(set_.G_) - 1)) ;
  }

  return df ;
}
